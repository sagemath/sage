/** @file numeric.h
 *
 *  Makes the interface to the underlying bignum package available. */

/*
 *  This is a modified version of code included with Ginac.  
 *  The modifications and modified version is:
 * 
 *      GiNaC-SAGE Copyright (C) 2008 William Stein
 *                           (C) 2015-16 Ralf Stephan
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


/*  The original copyright:
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

#ifndef __GINAC_NUMERIC_H__
#define __GINAC_NUMERIC_H__

#include "basic.h"
#include "constant.h"
#include "ex.h"

#include <gmp.h>
#include <limits>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <set>

void ginac_pyinit_Integer(PyObject*);
void ginac_pyinit_Float(PyObject*);
void ginac_pyinit_I(PyObject*);
PyObject* CC_get();

class CanonicalForm;

namespace GiNaC {

enum Type {
	LONG=1,
	PYOBJECT,
	MPZ,
	MPQ,
//	MPFR,
//	MPFC,
//	MPQC
};

union Value {
        Value() {}
        Value(signed long int i) : _long(i) {}
	signed long int _long;
	mpz_t _bigint;
	mpq_t _bigrat;
	PyObject* _pyobject;
};

extern const ex _ex0;

/** Exception class thrown when a singularity is encountered. */
class pole_error : public std::domain_error {
public:
	explicit pole_error(const std::string& what_arg, int degree);
	numeric degree() const;
private:
	int deg;
};

/** This class is a wrapper around CLN-numbers within the GiNaC class
 *  hierarchy. Objects of this type may directly be created by the user.*/
class numeric : public basic {
	GINAC_DECLARE_REGISTERED_CLASS(numeric, basic)
	// member functions

	// other constructors
public:
	numeric(const numeric&);
	numeric(int i) : basic(&numeric::tinfo_static), t(LONG), v(i)
        {
                hash = (i<0) ? i-1 : i;
                setflag(status_flags::evaluated | status_flags::expanded);
        }
	numeric(unsigned int i) : 
                basic(&numeric::tinfo_static), t(LONG), v(i), hash(i)
        { setflag(status_flags::evaluated | status_flags::expanded); }
	numeric(long i) : basic(&numeric::tinfo_static), t(LONG), v(i)
        {
                hash = (i<0) ? i-1 : i;
                setflag(status_flags::evaluated | status_flags::expanded);
        }
	numeric(unsigned long i) : 
                basic(&numeric::tinfo_static), t(LONG), v(i), hash(i)
        { setflag(status_flags::evaluated | status_flags::expanded); }
	numeric(long numer, long denom);
	numeric(double d);
	numeric(mpz_t bigint);
	numeric(mpq_t bigrat);
	numeric(PyObject*, bool=false);
        static ex unarchive(const archive_node &n, lst &sym_lst)
        {
                return (new numeric(n, sym_lst))->
                        setflag(status_flags::dynallocated);
        }

	~numeric();

	friend std::ostream& operator<<(std::ostream& os, const numeric& s);
	friend void coerce(numeric& new_left, numeric& new_right,
                        const numeric& left, const numeric& right);
        friend numeric & operator+=(numeric & lh, const numeric & rh);
        friend numeric & operator-=(numeric & lh, const numeric & rh);
        friend numeric & operator*=(numeric & lh, const numeric & rh);
        friend numeric & operator/=(numeric & lh, const numeric & rh);
        friend void rational_power_parts(const numeric& a, const numeric& b,
                        numeric& c, numeric& d, bool& c_unit);
public:

	unsigned precedence() const override { return 30; }
	bool info(unsigned inf) const override;
	bool is_polynomial(const ex & var) const override;
        bool is_long() const     { return t == LONG; }
	bool is_mpz() const      { return t == MPZ; }
	bool is_mpq() const      { return t == MPQ; }
        bool is_pyobject() const { return t == PYOBJECT; }
	bool is_zero() const;
	bool is_inexact_one() const;
	bool is_one() const;
	bool is_minus_one() const;
	bool is_positive() const override;
	bool is_negative() const override;
	bool is_integer() const override;
	bool is_pos_integer() const;
	bool is_nonneg_integer() const;
	bool is_even() const;
	bool is_odd() const;
	bool is_prime() const;
	bool is_rational() const;
	bool is_real() const override;
	bool is_cinteger() const;
	bool is_crational() const;
	bool is_exact() const;
        bool is_small_power(std::pair<int,int>& p) const;

	// functions overriding virtual functions from base classes
        numeric degree(const ex & s) const override;
	numeric ldegree(const ex & s) const override;
	ex coeff(const ex & s, const ex & n) const override;
	bool has(const ex &other, unsigned options = 0) const override;
	ex eval(int level = 0) const override;
	ex evalf(int level = 0, PyObject* parent = nullptr) const override;

	ex subs(const exmap & m, unsigned options = 0) const override;
	ex normal(exmap & repl, exmap & rev_lookup, int level = 0, unsigned options = 0) const override;
	ex to_rational(exmap & repl) const override;
	ex to_polynomial(exmap & repl) const override;
	numeric integer_content() const override;
	numeric max_coefficient() const override;
	numeric conj() const;
        ex conjugate() const override { return ex(conj()); }
	ex real_part() const override;
	ex imag_part() const override;
        ex series(const relational & r, int order, unsigned options) const override;

	/** Implementation of ex::diff for a numeric always returns 0.
	 *  @see ex::diff */
	ex derivative(const symbol &s) const override { return _ex0; }
        void useries(flint_series_t& fp, int order) const override;
	bool is_equal_same_type(const basic &other) const override;
	int compare_same_type(const numeric& right) const;
	long calchash() const override;

	const numeric add(const numeric &other) const;
	const numeric sub(const numeric &other) const;
	const numeric mul(const numeric &other) const;
	const numeric div(const numeric &other) const;
	const ex power(const numeric &other) const;
	const numeric power(signed long other) const;
        const numeric pow_intexp(const numeric &exponent) const;
	const numeric & add_dyn(const numeric &other) const;
	const numeric & sub_dyn(const numeric &other) const;
	const numeric & mul_dyn(const numeric &other) const;
	const numeric & div_dyn(const numeric &other) const;
	const numeric & power_dyn(const numeric &other) const;
	const numeric & operator=(int i);
	const numeric & operator=(unsigned int i);
	const numeric & operator=(long i);
	const numeric & operator=(double d);
	const numeric & operator=(const numeric& x);
        void swap(numeric& other);
	const numeric negative() const;
	const numeric inverse() const;
	const numeric step() const;
	int csgn() const;
	bool is_equal(const numeric &other) const;
	bool operator==(const numeric &other) const;
	bool operator!=(const numeric &other) const;
	bool operator<(const numeric &other) const;
	bool operator<=(const numeric &other) const;
	bool operator>(const numeric &other) const;
	bool operator>=(const numeric &other) const;
	int to_int() const;
	long to_long() const;
	double to_double() const;
        const numeric to_bigint() const;
        const mpz_t& as_mpz() const;
        const mpq_t& as_mpq() const;
        void canonicalize();
        PyObject* to_pyobject() const;
        const numeric try_py_method(const std::string& s) const;
        const numeric try_py_method(const std::string& s,
                        const numeric& x2) const;
        const numeric to_dict_parent(PyObject* dict) const;
        CanonicalForm to_canonical() const;

	const numeric real() const;
	const numeric imag() const;
	const numeric numer() const;
	const numeric denom() const;
        const numeric floor() const;
        const numeric frac() const;

        const numeric arbfunc_0arg(const char* name, PyObject* parent) const;
	const numeric exp(PyObject*) const;
	const numeric log(PyObject*) const;
	const numeric log(const numeric &b, PyObject*) const;
	const numeric ratlog(const numeric &b, bool& israt) const;
	const numeric sin() const;
	const numeric cos() const;
	const numeric tan() const;
	const numeric asin(PyObject*) const;
	const numeric acos(PyObject*) const;
	const numeric atan(PyObject*) const;
	const numeric atan(const numeric &y, PyObject*) const;
	const numeric sinh(PyObject*) const;
	const numeric cosh(PyObject*) const;
	const numeric tanh(PyObject*) const;
	const numeric asinh(PyObject*) const;
	const numeric acosh(PyObject*) const;
	const numeric atanh(PyObject*) const;
	const numeric Li2(const numeric &n, PyObject* parent) const;
	const numeric zeta() const;
	const numeric stieltjes() const;
	const numeric lgamma(PyObject* parent) const;
	const numeric gamma(PyObject* parent) const;
	const numeric rgamma(PyObject* parent) const;
	const numeric psi(PyObject* parent) const;
	const numeric psi(const numeric &n) const;
	const numeric factorial() const;
	const numeric doublefactorial() const;
	const numeric binomial(const numeric &k) const;
	const numeric bernoulli() const;
	const numeric hypergeometric_pFq(const std::vector<numeric>& a, const std::vector<numeric>& b, PyObject* parent) const;
	const numeric fibonacci() const;
	const numeric isqrt() const;
	const numeric sqrt() const;
        bool is_square() const;
        const ex sqrt_as_ex() const;
	const numeric abs() const;
	const numeric mod(const numeric &b) const;
	const numeric _smod(const numeric &b) const;
	ex smod(const numeric &b) const override;
	const numeric irem(const numeric &b) const;
	const numeric iquo(const numeric &b) const;
	const numeric iquo(const numeric &b, numeric &r) const;
	const numeric gcd(const numeric &b) const;
	const numeric lcm(const numeric &b) const;
        void factor(std::vector<std::pair<numeric, int>>& factors, long range=0) const;
        void factorsmall(std::vector<std::pair<long, int>>& factors, long range=0) const;
        void divisors(std::set<int>& divs) const;
	
protected:
	void print_numeric(const print_context & c, const char *par_open,
		const char *par_close, const char *imag_sym,
		const char *mul_sym, unsigned level, bool latex) const;
	void do_print(const print_context & c, unsigned level) const override;
	void do_print_latex(const print_latex & c, unsigned level) const;
	void do_print_tree(const print_tree & c, unsigned level) const override;
	void do_print_python_repr(const print_python_repr & c, unsigned level) const override;
        void dbgprint() const override;

        static bool integer_rational_power(numeric& res,
                const numeric& a, const numeric& b);

public:
        static const numeric binomial(unsigned long n, unsigned long k);

    Type t;
    Value v;
    long hash;
    bool is_hashable = true;
};


// global constants

extern numeric I;

// global functions

const numeric exp(const numeric &x, PyObject* parent=nullptr);
const numeric log(const numeric &x, PyObject* parent=nullptr);
const numeric log(const numeric &x, const numeric &b, PyObject* parent=nullptr);
const numeric sin(const numeric &x);
const numeric cos(const numeric &x);
const numeric tan(const numeric &x);
const numeric asin(const numeric &x, PyObject* parent=nullptr);
const numeric acos(const numeric &x, PyObject* parent=nullptr);
const numeric atan(const numeric &x, PyObject* parent=nullptr);
const numeric atan(const numeric &y, const numeric &x, PyObject* parent=nullptr);
const numeric sinh(const numeric &x, PyObject* parent=nullptr);
const numeric cosh(const numeric &x, PyObject* parent=nullptr);
const numeric tanh(const numeric &x, PyObject* parent=nullptr);
const numeric asinh(const numeric &x, PyObject* parent=nullptr);
const numeric acosh(const numeric &x, PyObject* parent=nullptr);
const numeric atanh(const numeric &x, PyObject* parent=nullptr);
const numeric Li2(const numeric &x, PyObject* parent=nullptr);
const numeric Li2(const numeric &x, const numeric &n, PyObject* parent=nullptr);
const numeric stieltjes(const numeric &x);
const numeric zeta(const numeric &x);
const numeric lgamma(const numeric &x, PyObject* parent=nullptr);
const numeric gamma(const numeric &x, PyObject* parent=nullptr);
const numeric rgamma(const numeric &x, PyObject* parent=nullptr);
const numeric psi(const numeric &x, PyObject* parent=nullptr);
const numeric psi(const numeric &n, const numeric &x);
const numeric beta(const numeric &x, const numeric &y, PyObject* parent=nullptr);
const numeric factorial(const numeric &n);
const numeric doublefactorial(const numeric &n);
const numeric binomial(const numeric &n, const numeric &k);
const numeric bernoulli(const numeric &n);
const numeric hypergeometric_pFq(const std::vector<numeric>& a, const std::vector<numeric>& b, const numeric &z, PyObject* parent);
const numeric fibonacci(const numeric &n);
const numeric isqrt(const numeric &x);
const numeric sqrt(const numeric &x);
const numeric abs(const numeric &x);
const numeric mod(const numeric &a, const numeric &b);
const numeric smod(const numeric &a, const numeric &b);
const numeric irem(const numeric &a, const numeric &b);
const numeric iquo(const numeric &a, const numeric &b);
const numeric iquo(const numeric &a, const numeric &b, numeric &r);
const numeric gcd(const numeric &a, const numeric &b);
const numeric lcm(const numeric &a, const numeric &b);

// wrapper functions around member functions

inline ex pow(const numeric &x, const numeric &y)
{
	return x.power(y);
}

inline const numeric inverse(const numeric &x)
{
	return x.inverse();
}

inline numeric step(const numeric &x)
{
	return x.step();
}

inline int csgn(const numeric &x)
{
	return x.csgn();
}

inline bool is_zero(const numeric &x)
{
	return x.is_zero();
}

inline bool is_positive(const numeric &x)
{
	return x.is_positive();
}

inline bool is_negative(const numeric &x)
{
	return x.is_negative();
}

inline bool is_integer(const numeric &x)
{
	return x.is_integer();
}

inline bool is_pos_integer(const numeric &x)
{
	return x.is_pos_integer();
}

inline bool is_nonneg_integer(const numeric &x)
{
	return x.is_nonneg_integer();
}

inline bool is_even(const numeric &x)
{
	return x.is_even();
}

inline bool is_odd(const numeric &x)
{
	return x.is_odd();
}

inline bool is_prime(const numeric &x)
{
	return x.is_prime();
}

inline bool is_rational(const numeric &x)
{
	return x.is_rational();
}

inline bool is_real(const numeric &x)
{
	return x.is_real();
}

inline bool is_cinteger(const numeric &x)
{
	return x.is_cinteger();
}

inline bool is_crational(const numeric &x)
{
	return x.is_crational();
}

inline int to_int(const numeric &x)
{
	return x.to_int();
}

inline long to_long(const numeric &x)
{
	return x.to_long();
}

inline double to_double(const numeric &x)
{
	return x.to_double();
}

inline const numeric real(const numeric &x)
{
	return x.real();
}

inline const numeric imag(const numeric &x)
{
	return x.imag();
}

inline const numeric numer(const numeric &x)
{
	return x.numer();
}

inline const numeric denom(const numeric &x)
{
	return x.denom();
}

// numeric evaluation functions for class constant objects:

ex ConstantEvalf(unsigned serial, PyObject* parent = nullptr);
ex UnsignedInfinityEvalf(unsigned serial, PyObject* parent = nullptr);
ex InfinityEvalf(unsigned serial, PyObject* parent = nullptr);
ex NegInfinityEvalf(unsigned serial, PyObject* parent = nullptr);

inline bool is_a_python_object(const ex & x)
{
        return (is_exactly_a<numeric>(x)
                and ex_to<numeric>(x).is_pyobject());
}

class conversion_error : public std::runtime_error {
    public:
        conversion_error() : std::runtime_error("") {}
};

} // namespace GiNaC

#endif // ndef __GINAC_NUMERIC_H__
