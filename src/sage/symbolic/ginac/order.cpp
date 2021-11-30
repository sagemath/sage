/** @file order.cpp
 *
 *  Definitions of order used for printing. */

/*
 *   Copyright (C) 2011 Burcin Erocal <burcin@erocal.org>
 *   Copyright (C) 2011 Jean-Pierre Flori <flori@enst.fr>
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

#include "order.h"
#include "registrar.h"
#include "utils.h"

#include <cmath>
#include <iostream>

namespace GiNaC {

double numeric_to_double(const numeric& exp)
{
	if (exp.is_real())
		return exp.to_double();
	
        return std::sqrt(std::pow(exp.real().to_double(), 2) + 
                         std::pow(exp.imag().to_double(), 2));
}


// Methods to get type id's. These cannot be const static because the
// id's are set by another static initializer.
  
const tinfo_t & print_order::function_id() const
{
	static tinfo_t id = find_tinfo_key("function");
	return id;
}

const tinfo_t & print_order::fderivative_id() const
{
	static tinfo_t id = find_tinfo_key("fderivative");
	return id;
}

const tinfo_t & print_order::power_id() const
{
	static tinfo_t id = find_tinfo_key("power");
	return id;
}

const tinfo_t & print_order::symbol_id() const
{
	static tinfo_t id = find_tinfo_key("symbol");
	return id;
}

const tinfo_t & print_order::mul_id() const
{
	static tinfo_t id = find_tinfo_key("mul");
	return id;
}

const tinfo_t & print_order::add_id() const
{
	static tinfo_t id = find_tinfo_key("add");
	return id;
}

const tinfo_t & print_order::numeric_id() const
{
	static tinfo_t id = find_tinfo_key("numeric");
	return id;
}

const tinfo_t & print_order::constant_id() const
{
	static tinfo_t id = find_tinfo_key("constant");
	return id;
}

const tinfo_t & print_order::wildcard_id() const
{
	static tinfo_t id = find_tinfo_key("wildcard");
	return id;
}

const tinfo_t & print_order::pseries_id() const
{
	static tinfo_t id = find_tinfo_key("pseries");
	return id;
}


/** What Sage does for printing:
 To print multivariate polynomials, SAGE uses "reversed" (bigger terms first)
 'degrevlex' order i.e. "reversed" graded reversed lexicographic order, the
 lexicographic order is defined below.  The variables are ordered according to
 their "creation" order i.e PR.<x,a,t> = PolynomialRing(QQ) gives x > a > t and
 x+t+a is printed x + a + t
*/

/** Returns true if lhex > rhex for 'degrevlex' of Sage
 i.e. graded reversed lexicographic order
*/
bool print_order::operator() (const ex &lhex, const ex &rhex) const 
{
	return compare(lhex, rhex) == 1;
}

int print_order::compare(const ex &lhex, const ex &rhex) const 
{
	return compare(*lhex.bp, *rhex.bp);
}

bool print_order_pair_mul::operator() (const expair &lhex, const expair &rhex) const
{
	int cmpval = print_order_mul().compare(lhex.rest, rhex.rest);
	if (cmpval != 0) {
		return cmpval == 1;
	}
	return compare_degrees(lhex, rhex);
}

bool print_order_pair::operator() (const expair &lhex, const expair &rhex) const
{
	// compare rests
	int cmpval = print_order().compare(lhex.rest, rhex.rest);
	if (cmpval != 0) {
		return cmpval == 1;
	}
	return compare_degrees(lhex, rhex);
}

bool print_order_pair::compare_degrees(const expair &lhex, const expair &rhex) const
{
	// compare coeffs which are numerics
	double lh_deg = numeric_to_double(ex_to<numeric>(lhex.coeff));
	double rh_deg = numeric_to_double(ex_to<numeric>(rhex.coeff));

	return lh_deg > rh_deg;	
}

/** Comparison functions:
 They should implement 'degrevlex' of Sage
 i.e. graded reversed lexicographic order
 The lexicographic order is the "natural" one on strings:
 a > b > ... > x > y > ... ?
*/

/** Return values for comparison functions :
 compare(a,b) should return :
 -1 if a < b
 0 if a == b
 1 if a > b
 as <=> in Perl and GiNaC internal functions
*/
int print_order::compare(const basic &lh, const basic &rh) const 
{
	const tinfo_t typeid_lh = lh.tinfo();
	const tinfo_t typeid_rh = rh.tinfo();

	if (typeid_rh==typeid_lh)
		if (typeid_rh == mul_id())
			return compare_same_type_mul(static_cast<const mul&>(lh),
						     static_cast<const mul&>(rh));
		else if (typeid_rh == add_id())
			return compare_same_type_add(
						     static_cast<const add&>(lh),
						     static_cast<const add&>(rh));
		else if (typeid_rh == symbol_id())
			return compare_same_type_symbol(static_cast<const symbol&>(lh),
							static_cast<const symbol&>(rh));
		else if (typeid_rh == power_id())
			return compare_same_type_power(static_cast<const power&>(lh),
						       static_cast<const power&>(rh));
		else if (typeid_rh == function_id())
			return compare_same_type_function(static_cast<const function&>(lh),
							  static_cast<const function&>(rh));
		else if (typeid_rh == fderivative_id())
			return compare_same_type_fderivative(static_cast<const fderivative&>(lh),
							     static_cast<const fderivative&>(rh));
		else
			// using GiNaC functions by default
			return lh.compare_same_type(rh);
		
	// at present numerics are combined into overall_coefficient
	// even when hold parameter is used
	else if (typeid_lh == numeric_id())
	 	//print numerics before anything else
		return 1;
	else if (typeid_rh == numeric_id())
	 	//print numerics before anything else
		return -1;
	else if (typeid_lh == wildcard_id())
	 	//print wildcards before anything else (but numerics)
		return 1;
	else if (typeid_rh == wildcard_id())
	 	//print wildcards before anything else (but numerics)
		return -1;
	else if (typeid_lh == constant_id())
	 	//print constants before anything else (but numerics, wildcards)
		return 1;
	else if (typeid_rh == constant_id())
	 	//print constants before anything else (but numerics, wildcards)
		return -1;
	else if (typeid_lh == fderivative_id())
		//print fderivatives after everything else
		return -1;
	else if (typeid_rh == fderivative_id())
		//print fderivatives after everything else
		return 1;
	else if (typeid_lh == function_id())
		//print functions before fderivatives, after anything else
		return -1;
	else if (typeid_rh == function_id())
		//print functions before fderivatives, after anything else
		return 1;
	else if (typeid_lh == mul_id()) {
		if (typeid_rh == power_id())
			return compare_mul_power(
					static_cast<const mul&>(lh),
					static_cast<const power&>(rh));
		if (typeid_rh == symbol_id())
			return compare_mul_symbol(
					static_cast<const mul&>(lh),
					static_cast<const symbol&>(rh));
		else if (typeid_rh == add_id())
			return -compare_add_mul(
					static_cast<const add&>(rh),
					static_cast<const mul&>(lh));
		else return generic_compare(typeid_lh, typeid_rh);
                }
        else if (typeid_lh == add_id()) {
		if (typeid_rh == power_id())
			return compare_add_power(
					static_cast<const add&>(lh),
					static_cast<const power&>(rh));
		if (typeid_rh == symbol_id())
			return compare_add_symbol(
					static_cast<const add&>(lh),
					static_cast<const symbol&>(rh));
		else if (typeid_rh == mul_id())
			return compare_add_mul(
					static_cast<const add&>(lh),
					static_cast<const mul&>(rh));
		else return generic_compare(typeid_lh, typeid_rh);
                }
	else if (typeid_lh == power_id()) {
		if (typeid_rh == mul_id())
			return -compare_mul_power(
					static_cast<const mul&>(rh),
					static_cast<const power&>(lh));
		if (typeid_rh == add_id())
			return -compare_add_power(
					static_cast<const add&>(rh),
					static_cast<const power&>(lh));
		else if (typeid_rh == symbol_id())
			return compare_power_symbol(
					static_cast<const power&>(lh),
					static_cast<const symbol&>(rh));
		else return generic_compare(typeid_lh, typeid_rh);
                }
	else if (typeid_lh == symbol_id()) {
		if (typeid_rh == mul_id())
			return -compare_mul_symbol(
					static_cast<const mul&>(rh),
					static_cast<const symbol&>(lh));
		if (typeid_rh == add_id())
			return -compare_add_symbol(
					static_cast<const add&>(rh),
					static_cast<const symbol&>(lh));
		else if (typeid_rh == power_id())
			return -compare_power_symbol(
					static_cast<const power&>(rh),
					static_cast<const symbol&>(lh));
		else return generic_compare(typeid_lh, typeid_rh);
                }
	else if (typeid_lh == pseries_id())
		//print pseries after everything else
		return -1;
	else if (typeid_rh == pseries_id())
		//print pseries after everything else
		return 1;
        return generic_compare(typeid_lh, typeid_rh);
}

// compare a mul and a symbol objects
// same behavior within mul and add objects:
// the total degree of the symbol is 1
// check total degree of mul then compare smallest item to symbol
int print_order::compare_mul_symbol(const mul &lh, const symbol &rh) const
{
	int cmpval;

	double tdeg;
	tdeg = lh.total_degree();

	if (tdeg != 1)
		return tdeg > 1 ? 1 : -1;
	
	const expair smallest_item = lh.get_sorted_seq().back();

	// compare bases
	cmpval = compare(*smallest_item.rest.bp, rh);
	if (cmpval != 0) {
		return cmpval;
	}

	// compare exponents
	cmpval = -compare(*smallest_item.coeff.bp, *_num1_p);
	if (cmpval != 0) {
		return cmpval;
	}

	if (lh.seq.size() == 1 and lh.overall_coeff.is_one())
		return 0;

	// there is something of total degree 0 in front of the mul object
	return 1;
}

// compare a mul and a pow objects
// same behavior within mul and add objects:
// first we compare total degrees
// if equal we compare the smallest basis in the sequence to the basis in other
// then their exponents
int print_order::compare_mul_power(const mul &lh, const power &rh) const
{
	int cmpval;

	double lh_deg = lh.total_degree();
	double rh_deg = 1;
	numeric rh_exp;
	if (is_exactly_a<numeric>(rh.exponent)) {
		rh_deg = numeric_to_double(ex_to<numeric>(rh.exponent));
	}
	if (rh_deg != lh_deg)
		return lh_deg < rh_deg ? -1 : 1;
	// same total degree

	// smallest item is at the end of the sorted sequence
	const expair smallest_item = lh.get_sorted_seq().back();
	
	// compare bases
	cmpval = compare(smallest_item.rest, rh.basis);
	if (cmpval != 0) {
		return cmpval;
	}

	// compare exponents
	cmpval = -compare(smallest_item.coeff, rh.exponent);
	if (cmpval != 0) {
		return cmpval;
	}

	if (lh.seq.size() == 1 and lh.overall_coeff.is_one())
		return 0;

	// there is something of total degree 0 in front of the mul object
	return 1;
}

// compare two mul objects
// same behavior within mul and add objects:
// first we compare total degrees
// if equal we compare the basis of the smallest items
// then their exponents
// and so on
int print_order::compare_same_type_mul(const mul &lh, const mul &rh) const
{
	int cmpval;
	
	// compare total degrees
	double lh_deg = lh.total_degree();
	double rh_deg = rh.total_degree();
	if (lh_deg != rh_deg)
		return lh_deg < rh_deg ? -1 : 1;
	
	// compare each item in lh to corresponding element in rh
	const epvector & sorted_seq1 = lh.get_sorted_seq();
	const epvector & sorted_seq2 = rh.get_sorted_seq();
	auto cit1 = sorted_seq1.rbegin();
	auto cit2 = sorted_seq2.rbegin();
	auto last1 = sorted_seq1.rend();
	auto last2 = sorted_seq2.rend();

	for (; (cit1!=last1)&&(cit2!=last2); ++cit1, ++cit2) {
		// compare bases
	  cmpval = compare(cit1->rest, cit2->rest);
		if (cmpval != 0) {
			return cmpval;
		}

		// compare exponents
		cmpval = -compare(cit1->coeff, cit2->coeff);
		if (cmpval != 0) {
			return cmpval;
		}
	}

	// compare sizes
	if (cit1 != last1) 
		return 1;
	if (cit2 != last2)
		return -1;

	// compare overall_coeff
	cmpval = compare(lh.overall_coeff, rh.overall_coeff);

	return cmpval;
}

// compare an add and a symbol objects
// same behavior within mul and add objects:
// the coefficient of the symbol is 1
int print_order::compare_add_symbol(const add &lh, const symbol &rh) const
{
	int cmpval;

	const expair biggest_item = lh.get_sorted_seq().front();

	// compare bases
	cmpval = print_order().compare(*biggest_item.rest.bp, rh);
	if (cmpval != 0) {
		return cmpval;
	}

	// compare coefficients
	cmpval = compare(*biggest_item.coeff.bp, *_num1_p);
	if (cmpval != 0) {
		return cmpval;
	}

	if (lh.seq.size() == 1 and lh.overall_coeff.is_zero())
		return 0;

	// there is something at the end of the add object
	return 1;
}

// compare an add and a mul objects
// same behavior within mul and add objects:
// the coefficient of the mul object is 1
int print_order::compare_add_mul(const add &lh,
				 const mul &rh) const
{
	int cmpval;
	const expair biggest_item = lh.get_sorted_seq().front();

	// compare bases
	cmpval = print_order().compare(*biggest_item.rest.bp, rh);
	if (cmpval != 0) {
		return cmpval;
	}

	// compare coefficients
	cmpval = compare(*biggest_item.coeff.bp, *_num1_p);
	if (cmpval != 0) {
		return cmpval;
	}

	if (lh.seq.size() == 1 and lh.overall_coeff.is_zero())
		return 0;

	// there is something at the end of the object
	return 1;
}

// compare an add and a pow objects
// same behavior within mul and add objects:
// the coefficient of the power object is 1
int print_order::compare_add_power(const add &lh,
		const power &rh) const
{
	int cmpval;
	const expair biggest_item = lh.get_sorted_seq().front();

	// compare bases
	cmpval = print_order().compare(*biggest_item.rest.bp, rh);
	if (cmpval != 0) {
		return cmpval;
	}

	// compare coefficients
	cmpval = compare(*biggest_item.coeff.bp, *_num1_p);
	if (cmpval != 0) {
		return cmpval;
	}

	if (lh.seq.size() == 1 and lh.overall_coeff.is_zero())
		return 0;

	// there is something at the end of the object
	return 1;
}

// compare two add objects
// same behavior within mul and add objects:
// first we compare the basis of the biggest items
// then their coefficients
// and so on
int print_order::compare_same_type_add(const add &lh, const add &rh) const
{
	int cmpval;

	const epvector & sorted_seq1 = lh.get_sorted_seq();
	const epvector & sorted_seq2 = rh.get_sorted_seq();
	auto cit1 = sorted_seq1.begin();
	auto cit2 = sorted_seq2.begin();
	auto last1 = sorted_seq1.end();
	auto last2 = sorted_seq2.end();

	for (; (cit1!=last1)&&(cit2!=last2); ++cit1, ++cit2) {
		// compare bases
		cmpval = print_order().compare(cit1->rest, cit2->rest);
		if (cmpval != 0) {
			return cmpval;
		}

		// compare coefficients
		cmpval = compare(cit1->coeff, cit2->coeff);
		if (cmpval != 0) {
			return cmpval;
		}
	}

	// compare sizes
	if (cit1 != last1) 
		return 1;
	if (cit2 != last2)
		return -1;

	// compare overall_coeff
	cmpval = compare(lh.overall_coeff, rh.overall_coeff);
	
	return cmpval;
}

// compare a power and a symbol objects
// does not behave the same way inside add or mul
// in mul object:
// first we compare bases
// then exponents
// in add object:
// first exponents
// then bases
int print_order_mul::compare_power_symbol
(const power &lh, const symbol &rh) const
{
	int cmpval;
	cmpval = compare(*lh.basis.bp, rh);
	if (cmpval != 0)
		  return cmpval;
	
	cmpval = compare(*lh.exponent.bp, *_num1_p);
	
	return cmpval;
}

int print_order::compare_power_symbol(const power &lh,
		const symbol &rh) const
{
	int cmpval;

	// We are in an add object
	if (is_exactly_a<numeric>(lh.exponent)) {
		double lh_deg = numeric_to_double(ex_to<numeric>(lh.exponent));
		if (lh_deg != 1)
			return lh_deg < 1 ? -1 : 1;
	}

	cmpval = compare(*lh.basis.bp, rh);

	return cmpval;
}

// compare two power objects
// does not behave the same way inside add or mul
// in mul object:
// fist we compare bases
// then exponents
// in add object:
// first exponents
// then bases
int print_order_mul::compare_same_type_power(const power &lh, const power &rh) const
{
	int cmpval;
	cmpval = compare(lh.basis, rh.basis);
	if (cmpval != 0)
		  return cmpval;
	cmpval = compare(lh.exponent, rh.exponent);
	
	return cmpval;
}


int print_order::compare_same_type_power(const power &lh, const power &rh) const
{
	int cmpval;

	double lh_deg = 1;
	double rh_deg = 1;
	if (is_exactly_a<numeric>(lh.exponent)) {
		lh_deg = numeric_to_double(ex_to<numeric>(lh.exponent));
	}
	if (is_exactly_a<numeric>(rh.exponent)) {
		rh_deg = numeric_to_double(ex_to<numeric>(rh.exponent));
	}
	if (rh_deg != lh_deg)
		return lh_deg < rh_deg ? -1 : 1;

	cmpval = compare(lh.basis, rh.basis);
	if (cmpval != 0)
		return cmpval;

	if (is_exactly_a<numeric>(lh.exponent) && is_exactly_a<numeric>(rh.exponent))
		return 0;
	return compare(lh.exponent, rh.exponent);
}

// compare two symbol objects
// same behavior within mul and add objects:
// we compare names
int print_order::compare_same_type_symbol(const symbol &lh, const symbol &rh) const
{
	/* Reversed ordering on char encoding (i.e. "a" < "b") then length
	 * (i.e. "abc" < "abcd")
	   i.e. "x" > "y" and "xyz" > "xyzt" */
	if (lh.serial==rh.serial) return 0;
	return lh.name < rh.name ? 1 : -1;
}

// compare two containers of the same type
template <template <class T, class = std::allocator<T> > class C>
int print_order::compare_same_type_container(const container<C> &lh,
					     const container<C> &rh) const
{
	auto it1 = lh.seq.begin(), it1end = lh.seq.end(),
			      it2 = rh.seq.begin(), it2end = rh.seq.end();

	while (it1 != it1end && it2 != it2end) {
		int cmpval = compare(*it1, *it2);
		if (cmpval)
			return cmpval;
		++it1; ++it2;
	}

	return (it1 == it1end) ? (it2 == it2end ? 0 : -1) : 1;
}

// compare two function objects
// same behavior within mul and add objects:
// we compare names
int print_order::compare_same_type_function(const function &lh,
					    const function &rh) const
{

	if (lh.serial==rh.serial)
		return compare_same_type_container(lh, rh);
	return lh.get_name() < rh.get_name() ? 1 : -1;	

}

// compare two fderivative objects
// same behavior within mul and add objects:
// we compare names
int print_order::compare_same_type_fderivative(const fderivative &lh,
					       const fderivative &rh) const
{
	int cmpval = compare_same_type_function(lh, rh);
	if (cmpval != 0)
		return cmpval;
	if (lh.parameter_set != rh.parameter_set)
		return lh.parameter_set < rh.parameter_set ? 1 : -1;
	return 0;
}

} // namespace GiNaC
