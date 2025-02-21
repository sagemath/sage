/** @file useries.cpp
 *
 *  Functions for extended truncated univariate power series. */

/*
 *  Copyright (C) 2016  Ralf Stephan <ralf@ark.in-berlin.de>
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

#ifndef __PYNAC_USERIES_CPP__
#define __PYNAC_USERIES_CPP__

#include "useries.h"
#include "useries-flint.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "symbol.h"
#include "function.h"
#include "relational.h"
#include "inifcns.h"
#include "utils.h"

#include <unordered_map>
#include <unordered_set>

namespace GiNaC {

// Normalize series if offset positive.
static void normalize(flint_series_t& fp)
{
        if (fp.offset > 0) {
                fmpq_poly_shift_left(fp.ft, fp.ft, fp.offset);
                fp.offset = 0;
        }
}

// Check that constant coeff of series is zero.
static void check_poly_ccoeff_zero(const flint_series_t& fp)
{
        if (fp.offset > 0)
                return;
        if (fp.offset < 0)
                throw flint_error();
        fmpq_t c;
        fmpq_init(c);
        fmpq_poly_get_coeff_fmpq(c, fp.ft, 0);
        if (not fmpq_is_zero(c))
                throw flint_error();
        fmpq_clear(c);
}

// Check that constant coeff of series is one.
static void check_poly_ccoeff_one(const flint_series_t& fp)
{
        if (fp.offset != 0)
                throw flint_error();
        fmpq_t c;
        fmpq_init(c);
        fmpq_poly_get_coeff_fmpq(c, fp.ft, 0);
        if (not fmpq_is_one(c))
                throw flint_error();
        fmpq_clear(c);
}

// Return low degree of polynomial
long fmpq_poly_ldegree(const fmpq_poly_t& fp)
{
        if (fmpq_poly_is_zero(fp))
                return 0;
        long len = fmpq_poly_length(fp);
        for (long n=0; n<=len; n++) {
                fmpq_t c;
                fmpq_init(c);
                fmpq_poly_get_coeff_fmpq(c, fp, (slong)n);
                if (not fmpq_is_zero(c)) {
                        fmpq_clear(c);
                        return n;
                }
                fmpq_clear(c);
        }
        return 0;
}

static void exp_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_exp_series(fp.ft, arg.ft, order);
}

static void log_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_one(arg);
        fmpq_poly_log_series(fp.ft, arg.ft, order);
}

static void sin_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_sin_series(fp.ft, arg.ft, order);
}

static void cos_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_cos_series(fp.ft, arg.ft, order);
}

static void tan_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_tan_series(fp.ft, arg.ft, order);
}

static void cot_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_tan_series(fp.ft, arg.ft, order);
        long ldeg = fmpq_poly_ldegree(fp.ft);
        fmpq_poly_shift_right(fp.ft, fp.ft, ldeg);
        fmpq_poly_inv_series(fp.ft, fp.ft, order-ldeg);
        fp.offset = -ldeg;
}

static void sec_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_cos_series(fp.ft, arg.ft, order);
        long ldeg = fmpq_poly_ldegree(fp.ft);
        fmpq_poly_shift_right(fp.ft, fp.ft, ldeg);
        fmpq_poly_inv_series(fp.ft, fp.ft, order-ldeg);
        fp.offset = -ldeg;
}

static void csc_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_sin_series(fp.ft, arg.ft, order);
        long ldeg = fmpq_poly_ldegree(fp.ft);
        fmpq_poly_shift_right(fp.ft, fp.ft, ldeg);
        fmpq_poly_inv_series(fp.ft, fp.ft, order-ldeg);
        fp.offset = -ldeg;
}

static void asin_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_asin_series(fp.ft, arg.ft, order);
}

static void atan_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_atan_series(fp.ft, arg.ft, order);
}

static void sinh_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_sinh_series(fp.ft, arg.ft, order);
}

static void cosh_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_cosh_series(fp.ft, arg.ft, order);
}

static void tanh_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_tanh_series(fp.ft, arg.ft, order);
}

static void coth_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_tanh_series(fp.ft, arg.ft, order);
        long ldeg = fmpq_poly_ldegree(fp.ft);
        fmpq_poly_shift_right(fp.ft, fp.ft, ldeg);
        fmpq_poly_inv_series(fp.ft, fp.ft, order-ldeg);
        fp.offset = -ldeg;
}

static void sech_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_cosh_series(fp.ft, arg.ft, order);
        long ldeg = fmpq_poly_ldegree(fp.ft);
        fmpq_poly_shift_right(fp.ft, fp.ft, ldeg);
        fmpq_poly_inv_series(fp.ft, fp.ft, order-ldeg);
        fp.offset = -ldeg;
}

static void csch_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_sinh_series(fp.ft, arg.ft, order);
        long ldeg = fmpq_poly_ldegree(fp.ft);
        fmpq_poly_shift_right(fp.ft, fp.ft, ldeg);
        fmpq_poly_inv_series(fp.ft, fp.ft, order-ldeg);
        fp.offset = -ldeg;
}

static void asinh_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_asinh_series(fp.ft, arg.ft, order);
}

static void atanh_useries(flint_series_t& fp, flint_series_t& arg, int order)
{
        check_poly_ccoeff_zero(arg);
        fmpq_poly_atanh_series(fp.ft, arg.ft, order);
}

using usfun_t = decltype(exp_useries);
using funcmap_t = std::unordered_map<unsigned int,usfun_t*>;

static funcmap_t& funcmap()
{
        static funcmap_t _funcmap = {{
                {exp_SERIAL::serial, &exp_useries},
                {log_SERIAL::serial, &log_useries},
                {sin_SERIAL::serial, &sin_useries},
                {cos_SERIAL::serial, &cos_useries},
                {tan_SERIAL::serial, &tan_useries},
                {cot_SERIAL::serial, &cot_useries},
                {sec_SERIAL::serial, &sec_useries},
                {csc_SERIAL::serial, &csc_useries},
                {asin_SERIAL::serial, &asin_useries},
                {atan_SERIAL::serial, &atan_useries},
                {sinh_SERIAL::serial, &sinh_useries},
                {cosh_SERIAL::serial, &cosh_useries},
                {tanh_SERIAL::serial, &tanh_useries},
                {coth_SERIAL::serial, &coth_useries},
                {sech_SERIAL::serial, &sech_useries},
                {csch_SERIAL::serial, &csch_useries},
                {asinh_SERIAL::serial, &asinh_useries},
                {atanh_SERIAL::serial, &atanh_useries},
        }};

        return _funcmap;
}

static bool rational_ex_f;

// Fast heuristic that rejects/accepts expressions for the fast
// expansion via Flint. It can give false positives that must be
// caught before Flint raises SIGABRT, because we want to use the
// older ::series() methods in case. Details:
//
// Does the expression have inexact values, constants, or such?
// It should practically consist of one symbol appearing in
// polynomials from QQ[], and only functions from a supported set.
// The helper uses recurrence to check that all numerics are from QQ,
// that there is not more than one symbol, no constants, and all
// function serial numbers are in the funcmap keys.
static bool unhandled_elements_in(const ex& the_ex, const symbol& symb)
{
        if (is_exactly_a<constant>(the_ex))
                return true;
        if (is_exactly_a<numeric>(the_ex))
                return not (ex_to<numeric>(the_ex).is_mpz()
                                or ex_to<numeric>(the_ex).is_long()
                                or ex_to<numeric>(the_ex).is_mpq());
        if (is_exactly_a<symbol>(the_ex)) {
                return (not ex_to<symbol>(the_ex).is_equal(symb));
        }
        if (is_exactly_a<function>(the_ex)) {
                rational_ex_f = false;
                const function& f = ex_to<function>(the_ex);
                if (funcmap().find(f.get_serial()) == funcmap().end())
                        return true;
                for (unsigned int i=0; i<f.nops(); i++)
                        if (unhandled_elements_in(f.op(i), symb))
                                return true;
                return false;
        }
        if (is_exactly_a<power>(the_ex)) {
                const power& pow = ex_to<power>(the_ex);
                if (not is_exactly_a<numeric>(pow.op(1)))
                        rational_ex_f = false;
                return (unhandled_elements_in(pow.op(0), symb)
                     or unhandled_elements_in(pow.op(1), symb));
        }
        if (is_a<expairseq>(the_ex)) {
                const expairseq& epseq = ex_to<expairseq>(the_ex);
                for (unsigned int i=0; i<epseq.nops(); i++) {
                        if (unhandled_elements_in(epseq.op(i), symb))
                                return true;
                }
                if (unhandled_elements_in(epseq.get_overall_coeff(), symb))
                        return true;
                return false;
        }
        return true;
}

bool useries_can_handle(const ex& the_ex, const symbol& s)
{
        rational_ex_f = true;
        bool ok = (not unhandled_elements_in(the_ex, s));
        if (ok) {
                ex nd = the_ex.numer_denom();
                try {
                        (void) nd.op(0).degree(s).to_long();
                        (void) nd.op(0).ldegree(s).to_long();
                        (void) nd.op(1).degree(s).to_long();
                        (void) nd.op(1).ldegree(s).to_long();
                }
                catch (conversion_error) {
                        throw std::runtime_error("exponent too big");
                }
                catch (std::runtime_error) {}
        }
        return ok;
}

class ldegree_error : public std::runtime_error {
    public:
        ldegree_error() : std::runtime_error("") {}
};

// This is practically ex::low_degree() except that some functions
// return nonzero values. We use this to determine if in the later
// series computation the precision may have to be increased. This
// is the case if we encounter an add in the treewalk. If not we
// can exactly determine the low degree.
static int low_series_degree(const ex& the_ex) {
        static std::unordered_set<unsigned int> funcset {{
                sin_SERIAL::serial,
                tan_SERIAL::serial,
                asin_SERIAL::serial,
                atan_SERIAL::serial,
                sinh_SERIAL::serial,
                tanh_SERIAL::serial,
                asinh_SERIAL::serial,
                atanh_SERIAL::serial,
}};

        if (is_exactly_a<constant>(the_ex)
            or is_exactly_a<numeric>(the_ex))
                return 0;
        if (is_exactly_a<symbol>(the_ex))
                return 1;
        if (is_exactly_a<function>(the_ex)) {
                const function& f = ex_to<function>(the_ex);
                unsigned int ser = f.get_serial();
                if (ser == log_SERIAL::serial)
                        return 1;
                if (ser == cot_SERIAL::serial
                    or ser == coth_SERIAL::serial
                    or ser == csc_SERIAL::serial
                    or ser == csch_SERIAL::serial)
                        return -low_series_degree(f.op(0));
                if (funcset.find(ser) == funcset.end())
                        return 0;
                return low_series_degree(f.op(0));
        }
        if (is_exactly_a<power>(the_ex)) {
                const power& pow = ex_to<power>(the_ex);
                ex expo = pow.op(1);
                if (is_exactly_a<numeric>(expo)) {
                        const numeric& n = ex_to<numeric>(expo);
                        if (n.is_integer())
                                return (low_series_degree(pow.op(0))
                                      * n.to_int());
                }
                return 0;
        }
        if (is_exactly_a<add>(the_ex)) {
                throw ldegree_error();
        }
        if (is_exactly_a<mul>(the_ex)) {
                int deg_sum = 0;
                const mul& m = ex_to<mul>(the_ex);
                for (const auto & elem : m.get_sorted_seq())
                        if (ex_to<numeric>(elem.coeff).is_integer())
                                deg_sum += low_series_degree(m.recombine_pair_to_ex(elem));
                return deg_sum;
        }
        return 0;
}

ex useries(const ex& the_ex, const symbol& x, int order, unsigned options)
{
        if (order <= 0)
                // send residues to the old code
                throw flint_error(); 
        bool may_extend = false;
        int ldeg = 0;
        try {
                ldeg = low_series_degree(the_ex);
        }
        catch (ldegree_error) {
                may_extend = true;
        }

        epvector epv;
        if (ldeg >= order) {
                epv.emplace_back(Order(_ex1), order);
                return pseries(relational(x,_ex0), epv);
        }

        if (ldeg > 0) {
                ldeg = 0;
                may_extend = false;
        }
        flint_series_t fp;
        fmpq_poly_set_ui(fp.ft, 0);
        int prec = order - ldeg + 2;
        the_ex.useries(fp, prec);
        int deg = fmpq_poly_degree(fp.ft);

        // Precision may have been lost when adding terms
        if (may_extend and deg < prec - fp.offset) {
                fmpq_poly_set_ui(fp.ft, 0);
                int old_offset = fp.offset;
                fp.offset = 0;
                the_ex.useries(fp, 2*prec - old_offset - deg);
        }

        // Fill expair vector
        for (int n=0; n<=deg+prec; n++) {
                if (n + fp.offset >= order)
                        break;
                fmpq_t c;
                fmpq_init(c);
                fmpq_poly_get_coeff_fmpq(c, fp.ft, (slong)n);
                if (not fmpq_is_zero(c)) {
                        mpq_t gc;
                        mpq_init(gc);
                        fmpq_get_mpq(gc, c);
                        numeric nc(gc); // numeric clears gc
                        epv.emplace_back(nc, numeric(n + fp.offset));
                }
                fmpq_clear(c);
        }
        epv.emplace_back(Order(_ex1), order);
        return pseries(relational(x,_ex0), epv);
}

void symbol::useries(flint_series_t& fp, int order) const
{
        fmpq_poly_set_str(fp.ft, "1  1");
        fp.offset = 1;
}

void add::useries(flint_series_t& fp, int order) const
{
        const numeric& oc = overall_coeff;
        if (oc.is_zero())
                fmpq_poly_set_ui(fp.ft, 0);
        else if (oc.is_long())
                fmpq_poly_set_si(fp.ft, oc.to_long());
        else if (oc.is_mpz()) {
                fmpz_t tmpfz;
                fmpz_init_set_readonly(tmpfz, oc.as_mpz());
                fmpq_poly_set_fmpz(fp.ft, tmpfz);
                fmpz_clear_readonly(tmpfz);
        } else {
                fmpq_t tmpfq;
                fmpq_init_set_readonly(tmpfq, oc.as_mpq());
                fmpq_poly_set_fmpq(fp.ft, tmpfq);
                fmpq_clear_readonly(tmpfq);
        }

        for (const auto & elem : seq) {
		const ex& t = recombine_pair_to_ex(elem);
                flint_series_t fp1;
                t.useries(fp1, order);
                if (fp.offset < fp1.offset) {
                        fmpq_poly_shift_left(fp1.ft, fp1.ft, fp1.offset-fp.offset);
                        fp1.offset = fp.offset;
                }
                else if (fp.offset > fp1.offset) {
                        fmpq_poly_shift_left(fp.ft, fp.ft, fp.offset-fp1.offset);
                        fp.offset = fp1.offset;
                }
                fmpq_poly_add(fp.ft, fp.ft, fp1.ft);
        }
}

void mul::useries(flint_series_t& fp, int order) const
{
        fmpq_poly_set_ui(fp.ft, 1);
        for (const auto & elem : seq) {
		const ex& t = recombine_pair_to_ex(elem);
                flint_series_t fp1;
                t.useries(fp1, order);
                long newoff = fp1.offset + fp.offset;
                fmpq_poly_mullow(fp.ft, fp.ft, fp1.ft, order+2);
                fp.offset = newoff;
        }
        const numeric& oc = overall_coeff;
        if (oc.is_one())
                return;

        if (oc.is_long())
                fmpq_poly_scalar_mul_si(fp.ft, fp.ft, oc.to_long());
        else if (oc.is_mpz()) {
                fmpz_t tmpfz;
                fmpz_init_set_readonly(tmpfz, oc.as_mpz());
                fmpq_poly_scalar_mul_fmpz(fp.ft, fp.ft, tmpfz);
                fmpz_clear_readonly(tmpfz);
        } else {
                fmpq_t tmpfq;
                fmpq_init_set_readonly(tmpfq, oc.as_mpq());
                fmpq_poly_scalar_mul_fmpq(fp.ft, fp.ft, tmpfq);
                fmpq_clear_readonly(tmpfq);
        }
}

void power::useries(flint_series_t& fp, int order) const
{
        flint_series_t fp1;
        basis.useries(fp1, order);
        if (not is_exactly_a<numeric>(exponent)) {
                check_poly_ccoeff_one(fp1);
                fmpq_poly_log_series(fp1.ft, fp1.ft, order);
                exponent.useries(fp, order);
                fmpq_poly_mullow(fp.ft, fp.ft, fp1.ft, order+2);
                check_poly_ccoeff_zero(fp);
                normalize(fp);
                fmpq_poly_exp_series(fp.ft, fp.ft, order);
                return;
        }
        const numeric& nexp = ex_to<numeric>(exponent);
        if (nexp.is_mpq()) {
                int num = nexp.numer().to_int();
                int den = nexp.denom().to_int();
                if (den == 2) { // exponent of form n/2
                        fmpq_t c;
                        fmpq_init(c);
                        fmpq_poly_get_coeff_fmpq(c, fp1.ft, 0);
                        mpz_t cnum, cden;
                        mpz_init(cnum);
                        mpz_init(cden);
                        fmpz_get_mpz(cnum, fmpq_numref(c));
                        fmpz_get_mpz(cden, fmpq_denref(c));
                        if (not mpz_perfect_square_p(cnum)
                            or not mpz_perfect_square_p(cden))
                                throw flint_error();
                        mpz_sqrt(cnum, cnum);
                        mpz_sqrt(cden, cden);
                        fmpq_t cc;
                        fmpz_init_set_readonly(fmpq_numref(cc), cnum);
                        fmpz_init_set_readonly(fmpq_denref(cc), cden);
                        mpz_clear(cnum);
                        mpz_clear(cden);

                        fmpq_poly_scalar_div_fmpq(fp1.ft, fp1.ft, c);
                        fmpq_poly_sqrt_series(fp1.ft, fp1.ft, order);
                        fmpq_poly_scalar_mul_fmpq(fp1.ft, fp1.ft, cc);
                        if (num > 0)
                                fmpq_poly_pow(fp.ft, fp1.ft, num);
                        else {
                                if (fmpq_poly_is_zero(fp1.ft))
                                        throw flint_error();
                                fmpq_poly_inv_series(fp1.ft, fp1.ft, order);
                                fmpq_poly_pow(fp.ft, fp1.ft, -num);
                        }
                        fmpq_clear(c);
                        return;
                }
                check_poly_ccoeff_one(fp1);
                fmpq_poly_log_series(fp1.ft, fp1.ft, order);
                fmpq_t tmp;
                fmpq_init_set_readonly(tmp, nexp.as_mpq());
                fmpq_poly_scalar_mul_fmpq(fp1.ft, fp1.ft, tmp);
                fmpq_clear_readonly(tmp);
                fmpq_poly_exp_series(fp.ft, fp1.ft, order);
                return;
        }
        // integer exponent
        int expint = nexp.to_int();
        long ldeg = fmpq_poly_ldegree(fp1.ft);
        if (expint > 0) {
                fmpq_poly_pow(fp.ft, fp1.ft, expint);
                fp.offset = fp1.offset * expint;
                fmpq_poly_truncate(fp.ft, fp.offset + order + 2);
                return;
        }
        if (expint < 0) {
                if (fmpq_poly_is_zero(fp1.ft))
                        throw flint_error();
                if (ldeg) {
                        fmpq_poly_shift_right(fp1.ft, fp1.ft, ldeg);
                        fp1.offset = ldeg;
                }
                fmpq_poly_inv_series(fp1.ft, fp1.ft,
                                order - fp1.offset*expint);
                fmpq_poly_pow(fp.ft, fp1.ft, -expint);
                fp.offset = fp1.offset * expint;
                fmpq_poly_truncate(fp.ft, order);
                return;
        }
        fmpq_poly_set_str(fp.ft, "1 1");
}

void function::useries(flint_series_t& fp, int order) const
{
        auto search = funcmap().find(serial);
        if (search == funcmap().end())
                throw std::runtime_error("can't happen in function::useries");
        flint_series_t fp1;
        seq[0].useries(fp1, order);
        normalize(fp1);
        (*search->second)(fp, fp1, order);
}

void numeric::useries(flint_series_t& fp, int order) const
{
        if (is_long())
                fmpq_poly_set_si(fp.ft, to_long());
        else if (is_mpz()) {
                fmpz_t tmpfz;
                fmpz_init_set_readonly(tmpfz, as_mpz());
                fmpq_poly_set_fmpz(fp.ft, tmpfz);
                fmpz_clear_readonly(tmpfz);
        } else {
                fmpq_t tmpfq;
                fmpq_init_set_readonly(tmpfq, as_mpq());
                fmpq_poly_set_fmpq(fp.ft, tmpfq);
                fmpq_clear_readonly(tmpfq);
        }
}

} // namespace GiNaC

#endif // undef __PYNAC_USERIES_CPP__
