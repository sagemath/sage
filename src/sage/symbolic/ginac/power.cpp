/** @file power.cpp
 *
 *  Implementation of GiNaC's symbolic exponentiation (basis^exponent). */

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

#include "power.h"
#include "expairseq.h"
#include "add.h"
#include "mul.h"
#include "numeric.h"
#include "constant.h"
#include "infinity.h"
#include "operators.h"
#include "inifcns.h" // for log() in power::derivative() and exp for printing
#include "symbol.h"
#include "lst.h"
#include "archive.h"
#include "utils.h"
#include "relational.h"
#include "compiler.h"
#include "cmatcher.h"
#include "wildcard.h"

#ifdef _WIN32
  #define NOMINMAX
  #include <windows.h>
#else
  #include <unistd.h>
#endif
#include <vector>
#include <stdexcept>
#include <limits>
#include <sstream>
#include <string>
#ifdef DO_GINAC_ASSERT
#  include <typeinfo>
#endif

namespace GiNaC {

GINAC_IMPLEMENT_REGISTERED_CLASS_OPT(power, basic,
  print_func<print_context>(&power::do_print).
  print_func<print_latex>(&power::do_print_latex).
  print_func<print_tree>(&basic::do_print_tree).
  print_func<print_python>(&power::do_print_python).
  print_func<print_python_repr>(&power::do_print_python_repr))

typedef std::vector<int> intvector;

//////////
// default constructor
//////////

power::power() : inherited(&power::tinfo_static) { }

//////////
// other constructors
//////////

// all inlined

//////////
// archiving
//////////

power::power(const archive_node &n, lst &sym_lst) : inherited(n, sym_lst)
{
	n.find_ex("basis", basis, sym_lst);
	n.find_ex("exponent", exponent, sym_lst);
}

void power::archive(archive_node &n) const
{
	inherited::archive(n);
	n.add_ex("basis", basis);
	n.add_ex("exponent", exponent);
}

//////////
// functions overriding virtual functions from base classes
//////////

// public

long power::calchash() const
{
  hashvalue = basis.gethash() + exponent.gethash();
  setflag(status_flags::hash_calculated);
  return hashvalue;
}

void power::print_power(const print_context & c, const char *powersymbol, const char *openbrace, const char *closebrace, unsigned level) const
{
	// Ordinary output of powers using '^' or '**'
	if (precedence() <= level) 
		c.s << openbrace << '(';
	basis.print(c, precedence());
	if (precedence() <= level)
		c.s << ')' << closebrace;
	c.s << powersymbol;
	c.s << openbrace;
	exponent.print(c, level);
	c.s << closebrace;
}

void power::do_print(const print_context & c, unsigned level) const
{
       if (exponent.is_minus_one()) {
		// inverses printed in a special way
 	        if (level >= 20) c.s << "(";
		c.s << "1/";
		basis.print(c, precedence());
 	        if (level >= 20) c.s << ")";
	} else if (exponent.is_equal(_ex1_2)) {

		// Square roots are printed in a special way
		c.s << "sqrt(";
		basis.print(c);
		c.s << ')';

	} else if (exponent.is_equal(_ex_1_2)) {
		// Square roots are printed in a special way
		c.s << "1/sqrt(";
		basis.print(c);
		c.s << ')';
	} else {

		std::stringstream tstream;
		print_context tcontext(tstream, c.options);
		exponent.print(tcontext, precedence());
		std::string expstr = tstream.str();
                if (expstr[0] == '-') {
			c.s<<"1/";
			expstr = expstr.erase(0, 1);
		}
		if (precedence() <= level)
			c.s << '(';
		// exp function prints as e^a. Printing powers of this can be
		// confusing, so we add parenthesis if the basis is exp
		bool exp_parenthesis = is_ex_the_function(basis, exp) and
                        not basis.op(0).is_one();
		if (exp_parenthesis)
			c.s << '(';
		basis.print(c, precedence());
		if (exp_parenthesis)
			c.s << ')';
		if (not exponent.is_minus_one()) {
    		        c.s << "^" << expstr;
		}
		if (precedence() <= level)
			c.s << ')';
	}
		
}

void power::do_print_latex(const print_latex & c, unsigned level) const
{
	if (is_exactly_a<numeric>(exponent) && ex_to<numeric>(exponent).is_negative()) {
		// Powers with negative numeric exponents are printed as fractions
		c.s << "\\frac{1}{";
		power(basis, -exponent).eval().print(c);
		c.s << '}';

	} else if (exponent.is_equal(_ex1_2)) {
		// Square roots are printed in a special way
		c.s << "\\sqrt{";
		basis.print(c);
		c.s << '}';

	} else {
		// exp function prints as e^a. Printing powers of this can be
		// confusing, so we add parenthesis if the basis is exp
		bool base_parenthesis = is_ex_the_function(basis, exp) and
			not basis.op(0).is_one();

		if (precedence() <= level)
			c.s << "{\\left(";
		if (base_parenthesis)
			c.s << "\\left(";

		basis.print(c, precedence());

		if (base_parenthesis)
			c.s << "\\right)";

		if (not exponent.is_minus_one()) {
    		        c.s << "^{";
			bool exp_parenthesis = is_exactly_a<power>(exponent);
			if (exp_parenthesis)
				c.s << "\\left(";
			exponent.print(c, 0);
			if (exp_parenthesis)
				c.s << "\\right)";
			c.s << '}';
		}

		if (precedence() <= level)
			c.s << "\\right)"<<'}';
	}
}

void power::do_print_python(const print_python & c, unsigned level) const
{
	print_power(c, "**", "", "", level);
}

void power::do_print_python_repr(const print_python_repr & c, unsigned level) const
{
	c.s << class_name() << '(';
	basis.print(c);
	c.s << ',';
	exponent.print(c);
	c.s << ')';
}

bool power::info(unsigned inf) const
{
        switch (inf) {
        case info_flags::polynomial:
        case info_flags::integer_polynomial:
        case info_flags::cinteger_polynomial:
        case info_flags::rational_polynomial:
        case info_flags::crational_polynomial:
        case info_flags::integer:
                return exponent.info(info_flags::nonnegint)
                and basis.info(inf);
        case info_flags::even:
                return exponent.info(info_flags::posint)
                and basis.is_integer();
        case info_flags::rational_function:
        case info_flags::rational:
                return exponent.is_integer()
                       and basis.info(inf);
        case info_flags::inexact:
                return exponent.info(inf)
                or basis.info(inf);
        case info_flags::algebraic:
                return !exponent.is_integer()
                or basis.info(inf);
        case info_flags::expanded:
                return (flags & status_flags::expanded) != 0u;
        case info_flags::positive:
                if (exponent.info(info_flags::even))
                        return basis.is_real()
                        and basis.info(info_flags::nonzero);
                if (exponent.info(info_flags::odd))
                        return basis.is_positive();
                return basis.is_positive()
                       and exponent.is_real();
        case info_flags::nonnegative:
                return (basis.is_positive()
                        and exponent.is_real())
                    or (basis.is_real()
                        and exponent.is_integer()
                        and exponent.info(info_flags::even));
        case info_flags::negative:
                if (exponent.info(info_flags::odd))
                        return basis.info(inf);
                return false;
        case info_flags::real:
                return ((basis.info(inf)
                         and exponent.is_integer())
                     or (basis.is_positive()
                         and exponent.info(inf)));
        case info_flags::nonzero:
                return (basis.info(inf)
                     or exponent.is_zero()
                     or exponent.info(info_flags::negative));
        }
        return inherited::info(inf);
}

size_t power::nops() const
{
	return 2;
}

const ex power::op(size_t i) const
{
	GINAC_ASSERT(i<2);

	return i==0 ? basis : exponent;
}

ex & power::let_op(size_t i)
{
	GINAC_ASSERT(i<2);

	return i==0 ? basis : exponent;
}

ex power::map(map_function & f) const
{
	const ex &mapped_basis = f(basis);
	const ex &mapped_exponent = f(exponent);

	if (!are_ex_trivially_equal(basis, mapped_basis)
	 || !are_ex_trivially_equal(exponent, mapped_exponent))
		return (new power(mapped_basis, mapped_exponent))->setflag(status_flags::dynallocated);
	
	return *this;
}

bool power::is_polynomial(const ex & var) const
{
	if (basis.is_polynomial(var)) {
		if (basis.has(var))
			// basis is non-constant polynomial in var
			return exponent.info(info_flags::nonnegint);
		
                // basis is constant in var
                return !exponent.has(var);
	}
	// basis is a non-polynomial function of var
	return false;
}

numeric power::degree(const ex & s) const
{
	if (is_equal(ex_to<basic>(s)))
		return *_num1_p;
	if (is_exactly_a<numeric>(exponent)
            and ex_to<numeric>(exponent).is_real()) {
		if (basis.is_equal(s))
			return ex_to<numeric>(exponent);
                return basis.degree(s) * ex_to<numeric>(exponent);
	} else if (basis.has(s))
		throw(std::runtime_error("power::degree(): undefined degree because of non-integer exponent"));
	else
		return *_num0_p;
}

numeric power::ldegree(const ex & s) const 
{
	if (is_equal(ex_to<basic>(s)))
		return *_num1_p;
	if (is_exactly_a<numeric>(exponent)
            and ex_to<numeric>(exponent).is_real()) {
		if (basis.is_equal(s))
			return ex_to<numeric>(exponent);
		return basis.ldegree(s) * ex_to<numeric>(exponent);
	} else if (basis.has(s))
		throw(std::runtime_error("power::ldegree(): undefined degree because of non-integer exponent"));
	else
		return *_num0_p;
}

ex power::coeff(const ex & s, const ex & n) const
{
	if (is_equal(ex_to<basic>(s)))
		return n.is_one() ? _ex1 : _ex0;
	if (!basis.is_equal(s)) {
		// basis not equal to s
		if (n.is_zero())
			return *this;
	} else {
		// basis equal to s
                if (not n.is_zero()
                    and exponent.is_equal(n))
                        return _ex1;
	}
        return _ex0;
}

/** Perform automatic term rewriting rules in this class.  In the following
 *  x, x1, x2,... stand for a symbolic variables of type ex and c, c1, c2...
 *  stand for such expressions that contain a plain number.
 *  - ^(x,0) -> 1  (also handles ^(0,0))
 *  - ^(x,1) -> x
 *  - ^(0,c) -> 0 or exception  (depending on the real part of c)
 *  - ^(1,x) -> 1
 *  - ^(c1,c2) -> *(c1^n,c1^(c2-n))  (so that 0<(c2-n)<1, try to evaluate roots, possibly in numerator and denominator of c1)
 *  - ^(^(x,c1),c2) -> ^(x,c1*c2)  if x is positive and c1 is real.
 *  - ^(^(x,c1),c2) -> ^(x,c1*c2)  (c2 integer or -1 < c1 <= 1 or (c1=-1 and c2>0), case c1=1 should not happen, see below!)
 *  - ^(*(x,y,z),c) -> *(x^c,y^c,z^c)  (if c integer)
 *  - ^(*(x,c1),c2) -> ^(x,c2)*c1^c2  (c1>0)
 *  - ^(*(x,c1),c2) -> ^(-x,c2)*c1^c2  (c1<0)
 *
 *  @param level cut-off in recursive evaluation */
ex power::eval(int level) const
{
	if (level == 1 and is_evaluated())
		return *this;
	if (level == -max_recursion_level)
		throw(std::runtime_error("max recursion level reached"));
	
	const ex & ebasis    = level==1 ? basis    : basis.eval(level-1);
	const ex & eexponent = level==1 ? exponent : exponent.eval(level-1);
	
	bool basis_is_numerical = false;
	bool exponent_is_numerical = false;
	numeric num_basis;
	numeric num_exponent;
	
	if (is_exactly_a<numeric>(ebasis)) {
		basis_is_numerical = true;
		num_basis = ex_to<numeric>(ebasis);

                // (negative num)^even --> (positive num)^even
                if (num_basis.is_negative()
                    and eexponent.info(info_flags::even)
                    and eexponent.is_real())
                        return power(-num_basis, eexponent);
	}
	if (is_exactly_a<numeric>(eexponent)) {
		exponent_is_numerical = true;
		num_exponent = ex_to<numeric>(eexponent);
	}

	// ^(\infty, x)
	if (is_exactly_a<infinity>(ebasis)) {
		const infinity & basis_inf = ex_to<infinity>(ebasis);
		if (eexponent.nsymbols()>0)
			throw(std::domain_error("power::eval(): pow(Infinity, f(x)) is not defined."));
		if (eexponent.is_zero())
			throw(std::domain_error("power::eval(): pow(Infinity, 0) is undefined."));
		if (eexponent.info(info_flags::negative))
			return _ex0;
		if (eexponent.is_positive()) {
			if (basis_inf.is_unsigned_infinity())
				return UnsignedInfinity;
			
				return mul(pow(basis_inf.get_direction(), eexponent), Infinity);
                }
		throw(std::domain_error("power::eval(): pow(Infinity, c)"
					" for constant of undetermined sign is not defined."));
	}

	// ^(x, \infty)
	if (is_exactly_a<infinity>(eexponent)) {
		const infinity & exp_inf = ex_to<infinity>(eexponent);
		if (exp_inf.is_unsigned_infinity())
			throw(std::domain_error("power::eval(): pow(x, unsigned_infinity) is not defined."));
		if (ebasis.nsymbols()>0) 
			throw(std::domain_error("power::eval(): pow(f(x), infinity) is not defined."));
		// x^(c*oo) --> (x^c)^(+oo)
		const ex abs_base = abs(pow(ebasis, exp_inf.get_direction()));
		if (abs_base > _ex1) {
			if (ebasis.is_positive())
				return Infinity;
			
				return UnsignedInfinity;
                        }
		if (abs_base < _ex1) return _ex0;
		if (abs_base == _ex1)
			throw(std::domain_error("power::eval(): pow(1, Infinity) is not defined."));
		throw(std::domain_error("power::eval(): pow(c, Infinity)"
					" for unknown magnitude |c| is not defined."));
	}
	
	// ^(x,0) -> 1  (0^0 also handled here)
	if (eexponent.is_zero())
		return _ex1;
	
	// ^(x,1) -> x
	if (eexponent.is_one())
		return ebasis;

	// ^(0,c1) -> 0 or exception  (depending on real value of c1)
	if (exponent_is_numerical and ebasis.is_zero() ) {
		if ((num_exponent.real()).is_zero())
			throw (std::domain_error("power::eval(): pow(0,I) is undefined"));
		else if ((num_exponent.real()).is_negative())
			throw (pole_error("power::eval(): division by zero",1));
		else
			return _ex0;
	}

	// ^(1,x) -> 1
	if (ebasis.is_one())
		return _ex1;

	// power of a function calculated by separate rules defined for this function
        if (is_exactly_a<function>(ebasis)) {
                const function& f = ex_to<function>(ebasis);
                if (f.get_serial() == exp_SERIAL::serial) {
                        if (is_exactly_a<numeric>(f.op(0))
                            or exponent_is_numerical)
                                return exp(mul(eexponent, f.op(0)));
                }
		return ex_to<function>(ebasis).power(eexponent);
        }

	// Turn (x^c)^d into x^(c*d) in the case that x is positive and c is real,
        // or if d is integer (and positive to preserve fraction output).
        // Secondly, if both c,d are numeric negative then cancel minuses.
	if (is_exactly_a<power>(ebasis)) {
                if (((eexponent.is_integer()
                      and eexponent.is_positive())
                    or (ebasis.op(0).is_positive()
                      and ebasis.op(1).is_real()))
                    or (ebasis.op(1).info(info_flags::odd)
                        and ebasis.op(0).is_real()
                        and (ebasis.op(1)*eexponent).is_integer()))
		        return power(ebasis.op(0), ebasis.op(1) * eexponent);
        }

        // (negative oc)^even --> (positive oc)^even
        if (is_exactly_a<mul>(ebasis)
            and ex_to<mul>(ebasis).overall_coeff.is_negative()
            and eexponent.info(info_flags::even)
            and eexponent.is_real())
                return power(-ebasis, eexponent);

	if (exponent_is_numerical) {
		// ^(c1,c2) -> c1^c2  (c1, c2 numeric(),
		// except if c1,c2 are rational, but c1^c2 is not)
                if (basis_is_numerical) {
                        const bool basis_is_crational = num_basis.is_crational();
                        const bool exponent_is_rational = num_exponent.is_rational();
                        const bool exponent_is_crational = num_exponent.is_crational();
                        if (not basis_is_crational
                            or not exponent_is_crational) {
                                // return a plain float
                                ex e = num_basis.power(num_exponent);
                                if (not is_exactly_a<numeric>(e))
                                        return (new power(ebasis, eexponent))->
                                                setflag(status_flags::dynallocated
                                                     | status_flags::evaluated);
                                return e;
                        }

                        if (exponent_is_rational)
                                return num_basis.power(num_exponent);
		        return this->hold();
                }

                // (2*x + 6*y)^(-4) -> 1/16*(x + 3*y)^(-4)
                if (is_exactly_a<add>(ebasis))
                        return ex_to<add>(ebasis).power(num_exponent);

                // ^(^(x,c1),c2) -> ^(x,c1*c2)
		// (c1, c2 numeric(), c2 integer or -1 < c1 <= 1 or
                // (c1=-1 and c2>0),
		// case c1==1 should not happen, see below!)
		if (is_exactly_a<power>(ebasis)) {
			const power & sub_power = ex_to<power>(ebasis);
			const ex & sub_basis = sub_power.basis;
			const ex & sub_exponent = sub_power.exponent;
			if (is_exactly_a<numeric>(sub_exponent)) {
				const numeric & num_sub_exponent = ex_to<numeric>(sub_exponent);
                                if (num_sub_exponent.is_negative()
                                    and num_exponent.is_negative())
                                        return power(power(sub_basis,
                                                           -num_sub_exponent),
                                                     -num_exponent);
				GINAC_ASSERT(num_sub_exponent!=numeric(1));
				if (num_exponent.is_integer() || (abs(num_sub_exponent) - (*_num1_p)).is_negative() 
						|| (num_sub_exponent == *_num_1_p && num_exponent.is_positive())) {
					return power(sub_basis,num_sub_exponent.mul(num_exponent));
				}
                                numeric pexp = num_sub_exponent * num_exponent;
                                if (pexp.is_integer()
                                    and num_sub_exponent.is_even()
                                    and sub_basis.is_real())
                                        return power(abs(sub_basis), pexp);
			}
		}
	
		// ^(*(...,x;c1),c2) -> *(^(*(...,x;1),c2),c1^c2)
                // (c1, c2 numeric(), c1>0)
		// ^(*(...,x;c1),c2) -> *(^(*(...,x;-1),c2),(-c1)^c2)
                // (c1, c2 numeric(), c1<0)
		if (is_exactly_a<mul>(ebasis)) {
                        if (num_exponent.is_integer())
                                // ^(*(x,y,z),c1) -> *(x^c1,y^c1,z^c1)
                                // (c1 integer)
                                return expand_mul(ex_to<mul>(ebasis),
                                                  num_exponent, 0);

			const mul & mulref = ex_to<mul>(ebasis);
			if (not mulref.overall_coeff.is_one()) {
				const numeric & num_coeff = mulref.overall_coeff;
				if (num_coeff.is_real()) {
					if (num_coeff.is_positive()) {
						auto mulp = new mul(mulref);
						mulp->overall_coeff = *_num1_p;
						mulp->setflag(status_flags::dynallocated);
						mulp->clearflag(status_flags::evaluated);
						mulp->clearflag(status_flags::hash_calculated);
						mulp->seq_sorted.resize(0);
						return (new mul(power(*mulp,exponent),
						                power(num_coeff, num_exponent)))->setflag(status_flags::dynallocated);
					} 
						GINAC_ASSERT(num_coeff.is_negative());
						if (not num_coeff.is_minus_one()) {
							auto mulp = new mul(mulref);
							mulp->overall_coeff = *_num_1_p;
							mulp->setflag(status_flags::dynallocated);
							mulp->clearflag(status_flags::evaluated);
							mulp->clearflag(status_flags::hash_calculated);
							mulp->seq_sorted.resize(0);
							return (new mul(power(*mulp, exponent),
							                power(abs(num_coeff), num_exponent)))->setflag(status_flags::dynallocated);
						}
					
				}
			}
                        else {
                                // (x*y^(m/n)*z)^(r/s) ---> y^t*(x*z)^(r/s), if t integer
                                exvector outer, inner;
                                for (size_t i=0; i<mulref.nops(); ++i) {
                                        const ex& fac = mulref.op(i);
                                        ex pfac = power(fac, eexponent);
                                        if (not is_exactly_a<power>(pfac)
                                            or ex_to<power>(pfac).exponent.info(info_flags::integer))
                                                outer.push_back(pfac);
                                        else
                                                inner.push_back(fac);
                                }
                                if (!outer.empty()) {
                                        ex outex, innex;
                                        if (outer.size() == 1)
                                                outex = outer[0];
                                        else
                                                outex = mul(outer).hold();
                                        if (!inner.empty()) {
                                                if (inner.size() == 1)
                                                        innex = inner[0];
                                                else
                                                        innex = mul(inner).hold();
                                                ex p = power(innex, eexponent).hold();
                                                return (new mul(outex, p))->setflag(status_flags::dynallocated | status_flags::evaluated);
                                        }
                                        return outex;
                                }
                        }
                }
	}

        if (is_exactly_a<mul>(eexponent)
            and ex_to<mul>(eexponent).get_overall_coeff().is_negative()) {
	        ex p = (new power(ebasis, -eexponent))->setflag(status_flags::dynallocated);
                if (ebasis.is_minus_one())
                        return p;
	        return (new power(p, _ex_1))->setflag(status_flags::dynallocated |
	                                               status_flags::evaluated);
        }


	// Reduce x^(c/log(x)) to exp(c) if x is positive
	if (ebasis.is_positive()) {
		if (eexponent.is_equal(1/log(ebasis)))
			return exp(log(basis)*exponent);
		if (is_exactly_a<mul>(eexponent) and
			 std::any_of(eexponent.begin(), eexponent.end(),
                                [ebasis](ex e)
                                { return e.is_equal(1/log(ebasis)); }))
					return exp(log(basis)*exponent);
	}
	
	if (are_ex_trivially_equal(ebasis,basis) &&
	    are_ex_trivially_equal(eexponent,exponent)) {
		return this->hold();
	}
	return (new power(ebasis, eexponent))->setflag(status_flags::dynallocated |
	                                               status_flags::evaluated);
}

ex power::evalf(int level, PyObject* parent) const
{
	ex ebasis;
	ex eexponent;
	
	if (level==1) {
		ebasis = basis;
		eexponent = exponent;
	} else if (level == -max_recursion_level) {
		throw(std::runtime_error("max recursion level reached"));
	} else {
		ebasis = basis.evalf(level-1, parent);
		if (!is_exactly_a<numeric>(exponent))
			eexponent = exponent.evalf(level-1, parent);
		else
			eexponent = exponent;
	}

	return power(ebasis,eexponent);
}

bool power::has(const ex & other, unsigned options) const
{
	if ((options & has_options::algebraic) == 0u)
		return basic::has(other, options);
	if (!is_a<power>(other))
		return basic::has(other, options);
	if (is_exactly_a<numeric>(exponent)
            and is_exactly_a<numeric>(other.op(1))) {
                if (!exponent.info(info_flags::integer)
                                || !other.op(1).info(info_flags::integer))
                        return basic::has(other, options);
                if (exponent.info(info_flags::posint)
                                && other.op(1).info(info_flags::posint)
                                && ex_to<numeric>(exponent)	> ex_to<numeric>(other.op(1))
                                && basis.match(other.op(0)))
                        return true;
                if (exponent.info(info_flags::negint)
                                && other.op(1).info(info_flags::negint)
                                && ex_to<numeric>(exponent) < ex_to<numeric>(other.op(1))
                                && basis.match(other.op(0)))
                        return true;
        }
	return basic::has(other, options);
}

bool power::match(const ex & pattern, exmap& map) const
{
	if (is_exactly_a<wildcard>(pattern)) {
                const auto& it = map.find(pattern);
                if (it != map.end())
		        return is_equal(ex_to<basic>(it->second));
		map[pattern] = *this;
		return true;
	} 
        if (not is_exactly_a<power>(pattern))
                return false;
        CMatcher cm(*this, pattern, map);
        const opt_exmap& m = cm.get();
        if (not m)
                return false;
        map = m.value();
        return true;
}

// from mul.cpp
extern bool tryfactsubs(const ex &, const ex &, int &, lst &);

ex power::subs(const exmap & m, unsigned options) const
{	
	const ex &subsed_basis = basis.subs(m, options);
	const ex &subsed_exponent = exponent.subs(m, options);

        if (!are_ex_trivially_equal(basis, subsed_basis)
	 || !are_ex_trivially_equal(exponent, subsed_exponent)) {
		ex p = power(subsed_basis, subsed_exponent);
		if (!is_exactly_a<power>(p)) {
		    // trac 30378 and 31530: do not over-substitute
	        return p;
	    }
                ex t = ex_to<power>(p).subs_one_level(m, options);
                if ((t-*this).is_zero())
                        return p;
                return t;
        }

	if ((options & subs_options::algebraic) == 0u)
		return subs_one_level(m, options);

        for (const auto & elem : m) {
		int nummatches = std::numeric_limits<int>::max();
		lst repls;
		if (tryfactsubs(*this, elem.first, nummatches, repls))
			return (ex_to<basic>((*this) * power(elem.second.subs(ex(repls),
                                subs_options::no_pattern) / elem.first.subs(ex(repls),
                                subs_options::no_pattern),
                                nummatches))).subs_one_level(m, options);
	}

	return subs_one_level(m, options);
}

ex power::conjugate() const
{
	// conjugate(pow(x,y))==pow(conjugate(x),conjugate(y)) unless on the
	// branch cut which runs along the negative real axis.
	if (basis.is_positive()) {
		ex newexponent = exponent.conjugate();
		if (are_ex_trivially_equal(exponent, newexponent)) {
			return *this;
		}
		return (new power(basis, newexponent))->setflag(status_flags::dynallocated);
	}
	if (exponent.is_integer()) {
		ex newbasis = basis.conjugate();
		if (are_ex_trivially_equal(basis, newbasis)) {
			return *this;
		}
		return (new power(newbasis, exponent))->setflag(status_flags::dynallocated);
	}
	return conjugate_function(*this).hold();
}

ex power::real_part() const
{
	if (is_exactly_a<numeric>(exponent)
	            and ((ex_to<numeric>(exponent).t == LONG) or (ex_to<numeric>(exponent).t == MPZ))) {
	        numeric nexp = ex_to<numeric>(exponent);
		ex basis_real = basis.real_part();
		if (basis_real == basis)
			return *this;
		symbol a("a"),b("b");
                a.set_domain(domain::real);
                b.set_domain(domain::real);
		ex result;
		if (nexp >= 0)
			result = power(a+I*b, nexp);
		else
			result = power(a/(a*a+b*b)-I*b/(a*a+b*b), -nexp);
		result = result.expand();
		result = result.real_part();
		result = result.subs(lst( a==basis_real, b==basis.imag_part() ));
		return result;
	}
	
	ex a = basis.real_part();
	ex b = basis.imag_part();
	ex c = exponent.real_part();
	ex d = exponent.imag_part();
	return power(abs(basis),c)*exp(-d*atan2(b,a))*cos(c*atan2(b,a)+d*log(abs(basis)));
}

ex power::imag_part() const
{
	if (is_exactly_a<numeric>(exponent)
	            and ((ex_to<numeric>(exponent).t == LONG) or (ex_to<numeric>(exponent).t == MPZ))) {
	        numeric nexp = ex_to<numeric>(exponent);
		ex basis_real = basis.real_part();
		if (basis_real == basis)
			return 0;
		symbol a("a"),b("b");
                a.set_domain(domain::real);
                b.set_domain(domain::real);
		ex result;
		if (nexp >= 0)
			result = power(a+I*b, nexp);
		else
			result = power(a/(a*a+b*b)-I*b/(a*a+b*b), -nexp);
		result = result.expand();
		result = result.imag_part();
		result = result.subs(lst( a==basis_real, b==basis.imag_part() ));
		return result;
	}
	
	ex a=basis.real_part();
	ex b=basis.imag_part();
	ex c=exponent.real_part();
	ex d=exponent.imag_part();
	return power(abs(basis),c)*exp(-d*atan2(b,a))*sin(c*atan2(b,a)+d*log(abs(basis)));
}

// protected

// protected

/** Implementation of ex::diff() for a power.
 *  @see ex::diff */
ex power::derivative(const symbol & s) const
{
	if (is_a<numeric>(exponent)) {
		// D(b^r) = r * b^(r-1) * D(b) (faster than the formula below)
		epvector newseq;
		newseq.reserve(2);
		newseq.emplace_back(basis, exponent - _ex1);
		newseq.emplace_back(basis.diff(s), _ex1);
		return mul(newseq, ex_to<numeric>(exponent));
	} 
	    // If the exponent is not a function of s, we have the following nice
	    // looking formula.  We use this to avoid getting ugly and hard to
	    // read output.
	    // D(b^e) = e * b^(e-1) * D(b)
	    ex ediff = exponent.diff(s);
	    if (ediff == 0) 
		return mul(mul(exponent, power(basis, exponent-1)), basis.diff(s));
	    
	    // The formula in the general case.
	    // D(b^e) = b^e * (D(e)*ln(b) + e*D(b)/b)
	    return mul(*this,
		       add(mul(ediff, log(basis)),
			   mul(mul(exponent, basis.diff(s)), power(basis, _ex_1))));
	
}

int power::compare_same_type(const basic & other) const
{
	GINAC_ASSERT(is_exactly_a<power>(other));
	const power &o = static_cast<const power &>(other);

	int cmpval = basis.compare(o.basis);
	if (cmpval != 0)
		return cmpval;
	
		return exponent.compare(o.exponent);
}

unsigned power::return_type() const
{
	return basis.return_type();
}

tinfo_t power::return_type_tinfo() const
{
	return basis.return_type_tinfo();
}

ex power::expand(unsigned options) const
{
	if (is_a<symbol>(basis)
            and exponent.is_integer()) {
		// A special case worth optimizing.
		setflag(status_flags::expanded);
		return *this;
	}

	// (x*p)^c -> x^c * p^c, if p>0
	// makes sense before expanding the basis
	if (is_exactly_a<mul>(basis)) {
		const mul &m = ex_to<mul>(basis);
		exvector prodseq;
		epvector powseq;
		prodseq.reserve(m.seq.size() + 1);
		powseq.reserve(m.seq.size() + 1);
		bool possign = true;

		// search for positive/negative factors
                for (const auto & elem : m.seq) {
			const ex& e = m.recombine_pair_to_ex(elem);
			if (e.is_positive())
				prodseq.push_back(pow(e, exponent).expand(options));
// 			we delete the following 'else if' clause because it can lead to
// 			an infinite loop (see sagemath :trac:`30688`)
// 			TODO: find a bug-free treatment of negative factors
// 			else if (e.info(info_flags::negative)) {
// 				prodseq.push_back(pow(-e, exponent).expand(options));
// 				possign = !possign;
// 			}
			else
				powseq.push_back(elem);
		}

		// take care on the numeric coefficient
		ex coef = (possign? _ex1 : _ex_1);
		if (m.overall_coeff.is_positive()
                                and not m.overall_coeff.is_one())
			prodseq.emplace_back(power(m.overall_coeff, exponent));
		else if (m.overall_coeff.is_negative()
                         and not m.overall_coeff.is_minus_one())
			prodseq.emplace_back(power(-m.overall_coeff, exponent));
		else
			coef *= m.overall_coeff;

		// If positive/negative factors are found, then extract them.
		// In either case we set a flag to avoid the second run on a part
		// which does not have positive/negative terms.
		if (!prodseq.empty()) {
			ex newbasis = coef*mul(powseq);
			ex_to<basic>(newbasis).setflag(status_flags::purely_indefinite);
			return ((new mul(prodseq))->setflag(status_flags::dynallocated)*(new power(newbasis, exponent))->setflag(status_flags::dynallocated).expand(options)).expand(options);
		} 
			ex_to<basic>(basis).setflag(status_flags::purely_indefinite);
	}

	const ex expanded_basis = basis.expand(options);
	const ex expanded_exponent = exponent.expand(options);
	
	// x^(a+b) -> x^a * x^b
	if (is_exactly_a<add>(expanded_exponent)) {
		const add &a = ex_to<add>(expanded_exponent);
		exvector distrseq;
		distrseq.reserve(a.seq.size() + 1);
                for (const auto & elem : a.seq)
			distrseq.emplace_back(power(expanded_basis, a.recombine_pair_to_ex(elem)));
		
		// Make sure that e.g. (x+y)^(2+a) expands the (x+y)^2 factor
		if (a.overall_coeff.is_integer()) {
			const numeric &num_exponent = a.overall_coeff;
			int int_exponent = num_exponent.to_int();
			if (int_exponent > 0 && is_exactly_a<add>(expanded_basis))
				distrseq.push_back(expand_add(ex_to<add>(expanded_basis), int_exponent, options));
			else
				distrseq.emplace_back(power(expanded_basis, a.overall_coeff));
		} else
			distrseq.emplace_back(power(expanded_basis, a.overall_coeff));
		
		// Make sure that e.g. (x+y)^(1+a) -> x*(x+y)^a + y*(x+y)^a
		ex r = (new mul(distrseq))->setflag(status_flags::dynallocated);
		return r.expand(options);
	}
	
	if (!is_exactly_a<numeric>(expanded_exponent) ||
		!ex_to<numeric>(expanded_exponent).is_integer()) {
		if (are_ex_trivially_equal(basis,expanded_basis) && are_ex_trivially_equal(exponent,expanded_exponent)) {
			return this->hold();
		} 
			return (new power(expanded_basis,expanded_exponent))->setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0));
		
	}
	
	// integer numeric exponent
	const numeric & num_exponent = ex_to<numeric>(expanded_exponent);
	int int_exponent = num_exponent.to_int();
	
	// (x+y)^n
	if (is_exactly_a<add>(expanded_basis)) {
                if (int_exponent == 1)
                        return expanded_basis;
                if ((options & expand_options::expand_only_numerators) == 0 and
                        int_exponent == -1)
                        return dynallocate<power>(expanded_basis, _ex_1).
                             setflag(status_flags::expanded|status_flags::evaluated);
                if ((options & expand_options::expand_only_numerators) != 0 and
                        int_exponent < 0)
                        return this->hold();
                if (int_exponent >= 0 or
                        (options & expand_options::expand_only_numerators) != 0)
		        return expand_add(ex_to<add>(expanded_basis),
                                        int_exponent, options);
                
                        return dynallocate<power>(expand_add(ex_to<add>(expanded_basis),
                                        -int_exponent, options), _ex_1).
                        setflag(status_flags::expanded|status_flags::evaluated);
        }
	
	// (x*y)^n -> x^n * y^n
	if (is_exactly_a<mul>(expanded_basis))
		return expand_mul(ex_to<mul>(expanded_basis), num_exponent, options, true);
	
	// cannot expand further
	if (are_ex_trivially_equal(basis,expanded_basis) && are_ex_trivially_equal(exponent,expanded_exponent))
		return this->hold();
	
		return (new power(expanded_basis,expanded_exponent))->setflag(status_flags::dynallocated | (options == 0 ? status_flags::expanded : 0));
}

//////////
// new virtual functions which can be overridden by derived classes
//////////

// none

//////////
// non-virtual functions in this class
//////////

namespace {  // anonymous namespace for power::expand_add() helpers

/** Helper class to generate all bounded combinatorial partitions of an integer
 *  n with exactly m parts (including zero parts) in non-decreasing order.
 */
class partition_generator {
private:
	// Partitions n into m parts, not including zero parts.
	// (Cf. OEIS sequence A008284; implementation adapted from Jörg Arndt's
	// FXT library)
	struct mpartition2
	{
		// partition: x[1] + x[2] + ... + x[m] = n and sentinel x[0] == 0
		std::vector<int> x;
		int n;   // n>0
		int m;   // 0<m<=n
		mpartition2(unsigned n_, unsigned m_)
		  : x(m_+1), n(n_), m(m_)
		{
			for (int k=1; k<m; ++k)
				x[k] = 1;
			x[m] = n - m + 1;
		}
		bool next_partition()
		{
			int u = x[m];  // last element
			int k = m;
			int s = u;
			while (--k) {
				s += x[k];
				if (x[k] + 2 <= u)
					break;
			}
			if (k==0)
				return false;  // current is last
			int f = x[k] + 1;
			while (k < m) {
				x[k] = f;
				s -= f;
				++k;
			}
			x[m] = s;
			return true;
		}
	} mpgen;
	int m;  // number of parts 0<m<=n
	mutable std::vector<int> partition;  // current partition
public:
	partition_generator(unsigned n_, unsigned m_)
	  : mpgen(n_, 1), m(m_), partition(m_)
	{ }
	// returns current partition in non-decreasing order, padded with zeros
	const std::vector<int>& current() const
	{
		for (int i = 0; i < m - mpgen.m; ++i)
			partition[i] = 0;  // pad with zeros

		for (int i = m - mpgen.m; i < m; ++i)
			partition[i] = mpgen.x[i - m + mpgen.m + 1];

		return partition;
	}
	bool next()
	{
		if (!mpgen.next_partition()) {
			if (mpgen.m == m || mpgen.m == mpgen.n)
				return false;  // current is last
			// increment number of parts
			mpgen = mpartition2(mpgen.n, mpgen.m + 1);
		}
		return true;
	}
};

/** Helper class to generate all compositions of a partition of an integer n,
 *  starting with the compositions which has non-decreasing order.
 */
class composition_generator {
private:
	// Generates all distinct permutations of a multiset.
	// (Based on Aaron Williams' algorithm 1 from "Loopless Generation of
	// Multiset Permutations using a Constant Number of Variables by Prefix
	// Shifts." <http://webhome.csc.uvic.ca/~haron/CoolMulti.pdf>)
	struct coolmulti {
		// element of singly linked list
		struct element {
			int value;
			element* next;
			element(int val, element* n)
			  : value(val), next(n) {}
			~element()
			{   // recurses down to the end of the singly linked list
				delete next;
			}
		};
		element *head, *i, *after_i;
		// NB: Partition must be sorted in non-decreasing order.
		explicit coolmulti(const std::vector<int>& partition)
		  : head(nullptr), i(nullptr), after_i(nullptr)
		{
                        if (partition.empty())
                                throw std::runtime_error("can't happn in coolmuli");
                        head = new element(partition[0], head);
                        i = head;
			for (unsigned n = 1; n < partition.size(); ++n) {
				head = new element(partition[n], head);
				if (n <= 1)
					i = head;
			}
			after_i = i->next;
		}
		~coolmulti()
		{   // deletes singly linked list
			delete head;
		}
		void next_permutation()
		{
			element *before_k;
			if (after_i->next != nullptr && i->value >= after_i->next->value)
				before_k = after_i;
			else
				before_k = i;
			element *k = before_k->next;
			before_k->next = k->next;
			k->next = head;
			if (k->value < head->value)
				i = k;
			after_i = i->next;
			head = k;
		}
		bool finished() const
		{
			return after_i->next == nullptr && after_i->value >= head->value;
		}
	} cmgen;
	bool atend;  // needed for simplifying iteration over permutations
	bool trivial;  // likewise, true if all elements are equal
	mutable std::vector<int> composition;  // current compositions
public:
	explicit composition_generator(const std::vector<int>& partition)
	  : cmgen(partition), atend(false), trivial(true), composition(partition.size())
	{
		for (unsigned i=1; i<partition.size(); ++i)
			trivial = trivial && (partition[0] == partition[i]);
	}
	const std::vector<int>& current() const
	{
		coolmulti::element* it = cmgen.head;
		size_t i = 0;
		while (it != nullptr) {
			composition[i] = it->value;
			it = it->next;
			++i;
		}
		return composition;
	}
	bool next()
	{
		// This ugly contortion is needed because the original coolmulti
		// algorithm requires code duplication of the payload procedure,
		// one before the loop and one inside it.
		if (trivial || atend)
			return false;
		cmgen.next_permutation();
		atend = cmgen.finished();
		return true;
	}
};

/** Helper function to compute the multinomial coefficient n!/(p1!*p2!*...*pk!)
 *  where n = p1+p2+...+pk, i.e. p is a partition of n.
 */
const numeric
multinomial_coefficient(const std::vector<int> & p)
{
	numeric n = 0, d = 1;
	for (auto & it : p) {
		n += numeric(it);
		d *= factorial(numeric(it));
	}
	return factorial(n) / d;
}

}  // anonymous namespace

/** expand a^n where a is an add and n is a positive integer.
 *  @see power::expand */
ex power::expand_add(const add & a, long n, unsigned options) const
{
	// The special case power(+(x,...y;x),2) can be optimized better.
	if (n==2)
		return expand_add_2(a, options);

	// method:
	//
	// Consider base as the sum of all symbolic terms and the overall numeric
	// coefficient and apply the binomial theorem:
	// S = power(+(x,...,z;c),n)
	//   = power(+(+(x,...,z;0);c),n)
	//   = sum(binomial(n,k)*power(+(x,...,z;0),k)*c^(n-k), k=1..n) + c^n
	// Then, apply the multinomial theorem to expand all power(+(x,...,z;0),k):
	// The multinomial theorem is computed by an outer loop over all
	// partitions of the exponent and an inner loop over all compositions of
	// that partition. This method makes the expansion a combinatorial
	// problem and allows us to directly construct the expanded sum and also
	// to re-use the multinomial coefficients (since they depend only on the
	// partition, not on the composition).
	// 
	// multinomial power(+(x,y,z;0),3) example:
	// partition : compositions                : multinomial coefficient
	// [0,0,3]   : [3,0,0],[0,3,0],[0,0,3]     : 3!/(3!*0!*0!) = 1
	// [0,1,2]   : [2,1,0],[1,2,0],[2,0,1],... : 3!/(2!*1!*0!) = 3
	// [1,1,1]   : [1,1,1]                     : 3!/(1!*1!*1!) = 6
	//  =>  (x + y + z)^3 =
	//        x^3 + y^3 + z^3
	//      + 3*x^2*y + 3*x*y^2 + 3*y^2*z + 3*y*z^2 + 3*x*z^2 + 3*x^2*z
	//      + 6*x*y*z
	//
	// multinomial power(+(x,y,z;0),4) example:
	// partition : compositions                : multinomial coefficient
	// [0,0,4]   : [4,0,0],[0,4,0],[0,0,4]     : 4!/(4!*0!*0!) = 1
	// [0,1,3]   : [3,1,0],[1,3,0],[3,0,1],... : 4!/(3!*1!*0!) = 4
	// [0,2,2]   : [2,2,0],[2,0,2],[0,2,2]     : 4!/(2!*2!*0!) = 6
	// [1,1,2]   : [2,1,1],[1,2,1],[1,1,2]     : 4!/(2!*1!*1!) = 12
	// (no [1,1,1,1] partition since it has too many parts)
	//  =>  (x + y + z)^4 =
	//        x^4 + y^4 + z^4
	//      + 4*x^3*y + 4*x*y^3 + 4*y^3*z + 4*y*z^3 + 4*x*z^3 + 4*x^3*z
	//      + 6*x^2*y^2 + 6*y^2*z^2 + 6*x^2*z^2
	//      + 12*x^2*y*z + 12*x*y^2*z + 12*x*y*z^2
	//
	// Summary:
	// r = 0
	// for k from 0 to n:
	//     f = c^(n-k)*binomial(n,k)
	//     for p in all partitions of n with m parts (including zero parts):
	//         h = f * multinomial coefficient of p
	//         for c in all compositions of p:
	//             t = 1
	//             for e in all elements of c:
	//                 t = t * a[e]^e
	//             r = r + h*t
	// return r

	epvector result;
	// The number of terms will be the number of combinatorial compositions,
	// i.e. the number of unordered arrangements of m nonnegative integers
	// which sum up to n.  It is frequently written as C_n(m) and directly
	// related with binomial coefficients: binomial(n+m-1,m-1).
        long anops = a.nops() - 1;
	long result_size = numeric::binomial(n + anops, anops).to_long();
	if (not a.overall_coeff.is_zero()) {
		// the result's overall_coeff is one of the terms
		--result_size;
	}
	result.reserve(result_size);

	// Iterate over all terms in binomial expansion of
	// S = power(+(x,...,z;c),n)
	//   = sum(binomial(n,k)*power(+(x,...,z;0),k)*c^(n-k), k=1..n) + c^n
	for (int k = 1; k <= n; ++k) {
		numeric binomial_coefficient;  // binomial(n,k)*c^(n-k)
		if (a.overall_coeff.is_zero()) {
			// degenerate case with zero overall_coeff:
			// apply multinomial theorem directly to power(+(x,...z;0),n)
			binomial_coefficient = *_num1_p;
			if (k < n) {
				continue;
			}
		} else {
                        if (n == k)
                                binomial_coefficient = *_num1_p;
                        else
                                binomial_coefficient = a.overall_coeff.pow_intexp(n-k) * numeric::binomial(n, k);
		}

		// Multinomial expansion of power(+(x,...,z;0),k)*c^(n-k):
		// Iterate over all partitions of k with exactly as many parts as
		// there are symbolic terms in the basis (including zero parts).
		partition_generator partitions(k, a.seq.size());
		do {
			const std::vector<int>& partition = partitions.current();
			// All monomials of this partition have the same number of terms and the same coefficient.
			const unsigned msize = std::count_if(partition.begin(), partition.end(), [](int i) { return i > 0; });
			const numeric coef = multinomial_coefficient(partition) * binomial_coefficient;

			// Iterate over all compositions of the current partition.
			composition_generator compositions(partition);
			do {
				const std::vector<int>& the_exponent = compositions.current();
				epvector monomial;
				monomial.reserve(msize);
				numeric factor = coef;
				for (unsigned i = 0; i < the_exponent.size(); ++i) {
					const ex & r = a.seq[i].rest;
					GINAC_ASSERT(!is_exactly_a<add>(r));
					GINAC_ASSERT(!is_exactly_a<power>(r) ||
						     !is_exactly_a<numeric>(ex_to<power>(r).exponent) ||
						     !ex_to<numeric>(ex_to<power>(r).exponent).is_pos_integer() ||
						     !is_exactly_a<add>(ex_to<power>(r).basis) ||
						     !is_exactly_a<mul>(ex_to<power>(r).basis) ||
						     !is_exactly_a<power>(ex_to<power>(r).basis));
					GINAC_ASSERT(is_exactly_a<numeric>(a.seq[i].coeff));
					const numeric & c = ex_to<numeric>(a.seq[i].coeff);
					if (the_exponent[i] == 0) {
						// optimize away
					} else if (the_exponent[i] == 1) {
						// optimized
						monomial.emplace_back(r, _ex1);
						if (not c.is_one())
							factor = factor.mul(c);
					} else { // general case exponent[i] > 1
						monomial.emplace_back(r, the_exponent[i]);
						if (not c.is_one())
							factor = factor.mul(c.pow_intexp(the_exponent[i]));
					}
				}
				result.emplace_back(mul(std::move(monomial)).expand(options), factor);
			} while (compositions.next());
		} while (partitions.next());
	}


	GINAC_ASSERT(result.size() == result_size);
	if (a.overall_coeff.is_zero()) {
		return dynallocate<add>(std::move(result)).setflag(status_flags::expanded).eval();
	} 
		return dynallocate<add>(std::move(result), a.overall_coeff.power(n)).setflag(status_flags::expanded).eval();
	
}


/** Special case of power::expand_add. Expands a^2 where a is an add.
 *  @see power::expand_add */
ex power::expand_add_2(const add & a, unsigned options) const
{
	epvector result;
	size_t result_size = (a.nops() * (a.nops()+1)) / 2;
	if (!a.overall_coeff.is_zero()) {
		// the result's overall_coeff is one of the terms
		--result_size;
	}
	result.reserve(result_size);

	auto last = a.seq.end();

	// power(+(x,...,z;c),2)=power(+(x,...,z;0),2)+2*c*+(x,...,z;0)+c*c
	// first part: ignore overall_coeff and expand other terms
	for (auto cit0=a.seq.begin(); cit0!=last; ++cit0) {
		const ex & r = cit0->rest;
		const ex & c = cit0->coeff;
		const numeric& nc = ex_to<numeric>(c);

		GINAC_ASSERT(!is_exactly_a<add>(r));
		GINAC_ASSERT(!is_exactly_a<power>(r) ||
		             !is_exactly_a<numeric>(ex_to<power>(r).exponent) ||
		             !ex_to<numeric>(ex_to<power>(r).exponent).is_pos_integer() ||
		             !is_exactly_a<add>(ex_to<power>(r).basis) ||
		             !is_exactly_a<mul>(ex_to<power>(r).basis) ||
		             !is_exactly_a<power>(ex_to<power>(r).basis));
		
		if (c.is_one()) {
			if (is_exactly_a<mul>(r)) {
				result.emplace_back(expand_mul(ex_to<mul>(r), *_num2_p, options, true),
				                        _ex1);
			} else {
				result.emplace_back(dynallocate<power>(r, _ex2),
				                        _ex1);
			}
		} else {
			if (is_exactly_a<mul>(r)) {
				result.emplace_back(expand_mul(ex_to<mul>(r), *_num2_p, options, true),
				                        nc.pow_intexp(*_num2_p));
			} else {
				result.emplace_back(dynallocate<power>(r, _ex2),
				                        nc.pow_intexp(*_num2_p));
			}
		}

		for (auto cit1=cit0+1; cit1!=last; ++cit1) {
			const ex & r1 = cit1->rest;
			const ex & c1 = cit1->coeff;
			result.emplace_back(mul(r,r1).expand(options),
			                        _num2_p->mul(ex_to<numeric>(c)).mul_dyn(ex_to<numeric>(c1)));
		}
	}
	
	// second part: add terms coming from overall_coeff (if != 0)
	if (!a.overall_coeff.is_zero()) {
		for (auto & i : a.seq)
			result.push_back(a.combine_pair_with_coeff_to_pair(i, a.overall_coeff.mul(*_num2_p)));
	}

	GINAC_ASSERT(result.size() == result_size);

	if (a.overall_coeff.is_zero()) {
		return dynallocate<add>(std::move(result)).setflag(status_flags::expanded);
	} 
		return dynallocate<add>(std::move(result), a.overall_coeff.pow_intexp(*_num2_p)).setflag(status_flags::expanded);
	
}

/** Expand factors of m in m^n where m is a mul and n is an integer.
 *  @see power::expand */
ex power::expand_mul(const mul & m, const numeric & n, unsigned options, bool from_expand) const
{
	GINAC_ASSERT(n.is_integer());

	if (n.is_zero()) {
		return _ex1;
	}

	epvector distrseq;
	distrseq.reserve(m.seq.size());
	bool need_reexpand = false;

        for (const auto & elem : m.seq) {
		expair p = m.combine_pair_with_coeff_to_pair(elem, n);
		if (from_expand && is_exactly_a<add>(elem.rest) && ex_to<numeric>(p.coeff).is_pos_integer()) {
			// this happens when e.g. (a+b)^(1/2) gets squared and
			// the resulting product needs to be reexpanded
			need_reexpand = true;
		}
		distrseq.push_back(p);
	}

	const mul & result = static_cast<const mul &>((new mul(distrseq, m.overall_coeff.pow_intexp(n)))->setflag(status_flags::dynallocated));
	if (need_reexpand)
		return ex(result).expand(options);
	if (from_expand)
		return result.setflag(status_flags::expanded);
	return result;
}

} // namespace GiNaC
