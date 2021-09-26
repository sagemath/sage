/** @file mpoly-giac.cpp
 * 
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *  Copyright (C) 2016  Ralf Stephan
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

#include "pynac-config.h"

#ifdef PYNAC_HAVE_LIBGIAC

#include <string>
#include <iostream>
#include <sstream>
#undef _POSIX_C_SOURCE
#undef _XOPEN_SOURCE

#include "upoly.h"
#include "basic.h"
#include "ex.h"
#include "add.h"
#include "constant.h"
#include "expairseq.h"
#include "mul.h"
#include "numeric.h"
#include "power.h"
#include "operators.h"
#include "pseries.h"
#include "relational.h"
#include "symbol.h"
#include "function.h"
#include "utils.h"

#include <giac/global.h>
#include <giac/gausspol.h>

namespace GiNaC {

static unsigned int the_dimension = 7;

static giac::context * context_ptr=nullptr;
static giac::gen giac_zero;
static giac::gen giac_one;


inline giac::polynome gen2pol(const giac::gen& g) {
        return giac::polynome(giac::monomial<giac::gen>(g, the_dimension));
}

inline giac::gen num2gen(const numeric& n) {
        auto gp = n.to_giacgen(context_ptr);
        if (gp != nullptr)
                return std::move(*gp);
        else {
                std::stringstream ss;
                ss << n;
                return giac::gen(std::string(ss.str()), context_ptr);
        }
}

static giac::polynome replace_with_symbol(const ex& e, ex_int_map& map, exvector& revmap)
{
        // Expression already replaced? Then return the assigned symbol
        auto it = map.find(e);
        if (it != map.end()) {
                const int index = it->second;
                giac::monomial<giac::gen> mon(giac_one, index, the_dimension);
                return giac::polynome(mon);
        }

        // Otherwise create new symbol and add to dict
        const int index = revmap.size() + 1;
        map.insert(std::make_pair(e, index));
        revmap.push_back(e);
        giac::monomial<giac::gen> mon(giac_one, index, the_dimension);
        return giac::polynome(mon);
}

const giac::polynome basic::to_polynome(ex_int_map& map, exvector& revmap)
{
        std::cerr << *this << std::endl;
        throw std::runtime_error("basic::to_polynome: can't happen");
}

const giac::polynome num_to_polynome(const numeric& num,
                ex_int_map& map, exvector& revmap)
{
        if (num.is_real()) {
                if (num.is_integer() or num.is_rational())
                        return gen2pol(num2gen(num));
                else
                        return replace_with_symbol(num, map, revmap);
        } else { // complex
                numeric re = num.real();
                numeric im = num.imag();
                giac::polynome re_p(the_dimension);
                giac::polynome im_p(the_dimension);
                if (re.is_integer() or re.is_rational())
                        re_p = gen2pol(num2gen(re));
                else
                        re_p = replace_with_symbol(re, map, revmap);
                if (im.is_integer() or im.is_rational())
                        im_p = gen2pol(num2gen(im));
                else
                        im_p = replace_with_symbol(im, map, revmap);
                giac::polynome r = re_p + im_p * replace_with_symbol(I, map, revmap);
                return r;
        }
}

// Convert to giac polynomial over QQ, filling replacement dicts
// TODO: special case numeric mpz_t, int instead of string interface
const giac::polynome ex::to_polynome(ex_int_map& map, exvector& revmap) const
{
        if (is_exactly_a<add>(*this))
        {
                const add& a = ex_to<add>(*this);
                giac::polynome p = gen2pol(giac_zero);
                for (const auto& termex : a.seq) {
                        p = p + a.recombine_pair_to_ex(termex).to_polynome(map, revmap);
                }
                p = p + num_to_polynome(a.overall_coeff, map, revmap);
                return p;
        }
        else if (is_exactly_a<numeric>(*this))
        {
                return num_to_polynome(ex_to<numeric>(*this), map, revmap);
        }
        else if (is_exactly_a<mul>(*this))
        {
                const mul& m = ex_to<mul>(*this);
                giac::polynome p = gen2pol(giac_one);
                for (const auto& termex : m.seq) {
                        p *= m.recombine_pair_to_ex(termex).to_polynome(map, revmap);
                }
                p *= num_to_polynome(m.overall_coeff, map, revmap);
                return p;
        }
        else if (is_exactly_a<power>(*this))
        {
                const power& pow = ex_to<power>(*this);
                if (is_exactly_a<numeric>(pow.exponent)) {
                        numeric expo = ex_to<numeric>(pow.exponent);
                        if (pow.exponent.info(info_flags::posint))
                                return std::move(giac::pow(pow.basis.to_polynome(map, revmap), expo.to_int()));
                }
                return replace_with_symbol(*this, map, revmap);
        }

        return replace_with_symbol(*this, map, revmap);
}

static ex gen2ex(const giac::gen& gen)
{
        // we need to handle giac types _INT_, _ZINT, and _CPLX
        switch (gen.type) {
        case giac::_INT_:
                return numeric(gen.val);
        case giac::_ZINT:
                mpz_t bigint;
                mpz_init_set(bigint, *(gen.ref_ZINTptr()));
                return numeric(bigint);
        case giac::_CPLX:
                return (gen2ex(gen.ref_CPLXptr()[0])
                                + I * gen2ex(gen.ref_CPLXptr()[1]));
        case giac::_FRAC:
                return (gen2ex(gen.ref_FRACptr()->num)
                                / gen2ex(gen.ref_FRACptr()->den));
        default:
                std::ostringstream os;
                os << "gen2ex: can't happen: " << int(gen.type) << std::flush;
                throw std::runtime_error(os.str());
        }
}

static ex polynome_to_ex(const giac::polynome& p, const exvector& revmap)
{
        ex e = _ex0;
        for (const auto& mon : p.coord) {
                ex prod = _ex1;
                for (unsigned int varno=0; varno<mon.index.size(); ++varno) {
                        if (mon.index[varno] != 0)
                                prod *= power(revmap[varno], mon.index[varno]);
                }
                e += gen2ex(mon.value) * prod;
        }
        return e;
}

// GCD of two exes which are in polynomial form
ex gcdpoly(const ex &a, const ex &b, ex *ca=nullptr, ex *cb=nullptr, bool check_args=true)
{
        //std::cerr << "gcd(" << a << "," << b << ") = ";
        if (context_ptr == nullptr) {
                context_ptr=new giac::context();
                giac_zero = giac::gen(std::string("0"), context_ptr);
                giac_one = giac::gen(std::string("1"), context_ptr);
        }

        if (a.is_zero())
                return b;
        if (b.is_zero())
                return a;

        // GCD of numerics
	if (is_exactly_a<numeric>(a) && is_exactly_a<numeric>(b)) {
		numeric g = gcd(ex_to<numeric>(a), ex_to<numeric>(b));
		if ((ca != nullptr) || (cb != nullptr)) {
			if (g.is_zero()) {
				if (ca != nullptr)
					*ca = _ex0;
				if (cb != nullptr)
					*cb = _ex0;
			} else {
				if (ca != nullptr)
					*ca = ex_to<numeric>(a) / g;
				if (cb != nullptr)
					*cb = ex_to<numeric>(b) / g;
			}
		}
		return g;
	}

	// Partially factored cases (to avoid expanding large expressions)
	if (is_exactly_a<mul>(a)) {
		if (is_exactly_a<mul>(b) && b.nops() > a.nops())
			goto factored_b;
factored_a:
		size_t num = a.nops();
		exvector g; g.reserve(num);
		exvector acc_ca; acc_ca.reserve(num);
		ex part_b = b;
		for (size_t i=0; i<num; i++) {
			ex part_ca, part_cb;
			g.push_back(gcdpoly(a.op(i), part_b, &part_ca, &part_cb, check_args));
			acc_ca.push_back(part_ca);
			part_b = part_cb;
		}
		if (ca != nullptr)
			*ca = (new mul(acc_ca))->setflag(status_flags::dynallocated);
		if (cb != nullptr)
			*cb = part_b;
		return (new mul(g))->setflag(status_flags::dynallocated);
	} else if (is_exactly_a<mul>(b)) {
		if (is_exactly_a<mul>(a) && a.nops() > b.nops())
			goto factored_a;
factored_b:
		size_t num = b.nops();
		exvector g; g.reserve(num);
		exvector acc_cb; acc_cb.reserve(num);
		ex part_a = a;
		for (size_t i=0; i<num; i++) {
			ex part_ca, part_cb;
			g.push_back(gcdpoly(part_a, b.op(i), &part_ca, &part_cb, check_args));
			acc_cb.push_back(part_cb);
			part_a = part_ca;
		}
		if (ca != nullptr)
			*ca = part_a;
		if (cb != nullptr)
			*cb = (new mul(acc_cb))->setflag(status_flags::dynallocated);
		return (new mul(g))->setflag(status_flags::dynallocated);
	}

	// Input polynomials of the form poly^n are sometimes also trivial
	if (is_exactly_a<power>(a)) {
		ex p = a.op(0);
		const ex& exp_a = a.op(1);
		if (is_exactly_a<power>(b)) {
			ex pb = b.op(0);
			const ex& exp_b = b.op(1);
			if (p.is_equal(pb)) {
				// a = p^n, b = p^m, gcd = p^min(n, m)
				if (exp_a < exp_b) {
					if (ca != nullptr)
						*ca = _ex1;
					if (cb != nullptr)
						*cb = power(p, exp_b - exp_a);
					return power(p, exp_a);
				} else {
					if (ca != nullptr)
						*ca = power(p, exp_a - exp_b);
					if (cb != nullptr)
						*cb = _ex1;
					return power(p, exp_b);
				}
			} else {
				ex p_co, pb_co;
				ex p_gcd = gcdpoly(p, pb, &p_co, &pb_co, check_args);
				if (p_gcd.is_one()) {
					// a(x) = p(x)^n, b(x) = p_b(x)^m, gcd (p, p_b) = 1 ==>
					// gcd(a,b) = 1
					if (ca != nullptr)
						*ca = a;
					if (cb != nullptr)
						*cb = b;
					return _ex1;
					// XXX: do I need to check for p_gcd = -1?
				} else {
					// there are common factors:
					// a(x) = g(x)^n A(x)^n, b(x) = g(x)^m B(x)^m ==>
					// gcd(a, b) = g(x)^n gcd(A(x)^n, g(x)^(n-m) B(x)^m
					if (exp_a < exp_b) {
						return power(p_gcd, exp_a)*
							gcdpoly(power(p_co, exp_a), power(p_gcd, exp_b-exp_a)*power(pb_co, exp_b), ca, cb, false);
					} else {
						return power(p_gcd, exp_b)*
							gcdpoly(power(p_gcd, exp_a - exp_b)*power(p_co, exp_a), power(pb_co, exp_b), ca, cb, false);
					}
				} // p_gcd.is_equal(_ex1)
			} // p.is_equal(pb)

		} else {
			if (p.is_equal(b)) {
				// a = p^n, b = p, gcd = p
				if (ca != nullptr)
					*ca = power(p, a.op(1) - 1);
				if (cb != nullptr)
					*cb = _ex1;
				return p;
			} 

			ex p_co, bpart_co;
			ex p_gcd = gcdpoly(p, b, &p_co, &bpart_co, false);

			if (p_gcd.is_one()) {
				// a(x) = p(x)^n, gcd(p, b) = 1 ==> gcd(a, b) = 1
				if (ca != nullptr)
					*ca = a;
				if (cb != nullptr)
					*cb = b;
				return _ex1;
			} else {
				// a(x) = g(x)^n A(x)^n, b(x) = g(x) B(x) ==> gcd(a, b) = g(x) gcd(g(x)^(n-1) A(x)^n, B(x))
				return p_gcd*gcdpoly(power(p_gcd, exp_a-1)*power(p_co, exp_a), bpart_co, ca, cb, false);
			}
		} // is_exactly_a<power>(b)

	} else if (is_exactly_a<power>(b)) {
		ex p = b.op(0);
		if (p.is_equal(a)) {
			// a = p, b = p^n, gcd = p
			if (ca != nullptr)
				*ca = _ex1;
			if (cb != nullptr)
				*cb = power(p, b.op(1) - 1);
			return p;
		}

		ex p_co, apart_co;
		const ex& exp_b(b.op(1));
		ex p_gcd = gcdpoly(a, p, &apart_co, &p_co, false);
		if (p_gcd.is_one()) {
			// b=p(x)^n, gcd(a, p) = 1 ==> gcd(a, b) == 1
			if (ca != nullptr)
				*ca = a;
			if (cb != nullptr)
				*cb = b;
			return _ex1;
		} else {
			// there are common factors:
			// a(x) = g(x) A(x), b(x) = g(x)^n B(x)^n ==> gcd = g(x) gcd(g(x)^(n-1) A(x)^n, B(x))

			return p_gcd*gcdpoly(apart_co, power(p_gcd, exp_b-1)*power(p_co, exp_b), ca, cb, false);
		} // p_gcd.is_equal(_ex1)
	}

        // Conversion necessary to count needed symbols beforehand
        exmap repl;
        ex poly_a = a.to_rational(repl);
        ex poly_b = b.to_rational(repl);

        symbolset s1 = poly_a.symbols();
        const symbolset& s2 = poly_b.symbols();
        s1.insert(s2.begin(), s2.end());
        the_dimension = s1.size();

        ex_int_map map;
        exvector revmap;

        giac::polynome p = poly_a.to_polynome(map, revmap);
        giac::polynome q = poly_b.to_polynome(map, revmap);
        giac::polynome d(the_dimension);
        giac::gcd(p, q, d);

        if (ca != nullptr) {
                giac::polynome quo(the_dimension);
                if (giac::exactquotient(p, d, quo))
                        *ca = polynome_to_ex(quo, revmap).subs(repl, subs_options::no_pattern);
                else
                        throw(std::runtime_error("can't happen in gcdpoly"));
        }
        if (cb != nullptr) {
                giac::polynome quo(the_dimension);
                if (giac::exactquotient(q, d, quo))
                        *cb = polynome_to_ex(quo, revmap).subs(repl, subs_options::no_pattern);
                else
                        throw(std::runtime_error("can't happen in gcdpoly"));
        }
        return polynome_to_ex(d, revmap).subs(repl, subs_options::no_pattern);
}

#if 0
ex resultantpoly(const ex & ee1, const ex & ee2, const ex & s)
{
        // Conversion necessary to count needed symbols beforehand
        exmap repl;
        ex poly_a = ee1.to_rational(repl);
        ex poly_b = ee2.to_rational(repl);

        symbolset s1 = poly_a.symbols();
        const symbolset& s2 = poly_b.symbols();
        s1.insert(s2.begin(), s2.end());
        s1.insert(ex_to<symbol>(s));
        the_dimension = s1.size();

        ex_int_map map;
        exvector revmap;

        giac::polynome p = poly_a.to_polynome(map, revmap);
        giac::polynome q = poly_b.to_polynome(map, revmap);
        giac::polynome d = giac::resultant(p, q);
        return polynome_to_ex(d, revmap).subs(repl, subs_options::no_pattern);
}
#endif

bool factorpoly(const ex& the_ex, ex& res_prod)
{
        if (is_exactly_a<numeric>(the_ex)
            or is_exactly_a<function>(the_ex)
            or is_exactly_a<constant>(the_ex)
            or is_exactly_a<symbol>(the_ex))
                return false;

        if (is_exactly_a<mul>(the_ex)) {
                const mul& m = ex_to<mul>(the_ex);
                exvector ev;
                bool mchanged = false;
                for (size_t i=0; i<m.nops(); ++i) {
                        ex r;
                        const ex& e = m.op(i);
                        bool res = factorpoly(e, r);
                        if (res) {
                                ev.push_back(r);
                                mchanged = true;
                        }
                        else
                                ev.push_back(e);
                }
                if (mchanged)
                        res_prod = mul(ev);
                return mchanged;
        }

        if (is_exactly_a<power>(the_ex)) {
                const power& pow = ex_to<power>(the_ex);
                ex prod;
                bool res = factorpoly(pow.op(0), prod);
                if (not res)
                        return false;
                res_prod = power(prod, pow.op(1));
                return true;
        }

        if (not is_exactly_a<add>(the_ex))
                throw(std::runtime_error("can't happen in factor"));

        if (context_ptr == nullptr) {
                context_ptr=new giac::context();
                giac_zero = giac::gen(std::string("0"), context_ptr);
                giac_one = giac::gen(std::string("1"), context_ptr);
        }

        exmap repl;
        ex poly = the_ex.to_rational(repl);
        symbolset s1 = poly.symbols();
        the_dimension = s1.size();

        ex_int_map map;
        exvector revmap;
        giac::polynome p = the_ex.to_polynome(map, revmap);
        giac::polynome p_content(the_dimension);
        giac::factorization f;
        bool res = factor(p, p_content, f,
                        false, false, false, giac_one, giac_one);
        if (not res)
                return false;
        res_prod = polynome_to_ex(p_content, revmap).subs(repl, subs_options::no_pattern);
        for (auto fpair : f)
                res_prod = mul(res_prod,
                                power(polynome_to_ex(fpair.fact, revmap).subs(repl, subs_options::no_pattern),
                                        fpair.mult));
        return true;
}

ex poly_mul_expand(const ex& a, const ex& b)
{
        exmap repl;
        ex poly_a = a.to_rational(repl);
        ex poly_b = b.to_rational(repl);

        symbolset s1 = poly_a.symbols();
        const symbolset& s2 = poly_b.symbols();
        s1.insert(s2.begin(), s2.end());
        the_dimension = s1.size();

        ex_int_map map;
        exvector revmap;

        giac::polynome p = poly_a.to_polynome(map, revmap);
        giac::polynome q = poly_b.to_polynome(map, revmap);
        giac::polynome d(the_dimension);
        d = p * q;

        return polynome_to_ex(d, revmap).subs(repl, subs_options::no_pattern);
}

} // namespace GiNaC

#endif // HAVE_LIBGIAC
