/** @file inifcns_orthopoly.cpp
 *
 *  Implementation of orthogonal polynomials as functions.
 *
 *  (C) 2016 Ralf Stephan <ralf@ark.in-berlin.de>
 */

#include "inifcns.h"
#include "ex.h"
#include "constant.h"
#include "infinity.h"
#include "numeric.h"
#include "add.h"
#include "mul.h"
#include "power.h"
#include "operators.h"
#include "relational.h"
#include "symbol.h"
#include "pseries.h"
#include "utils.h"

#include "gmp.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpq_poly.h"
#include "flint/fmpq.h"

namespace GiNaC {

//////////
// Chebyshev polynomials T_n(x)
//////////

static ex chebyt_eval(const ex& n_, const ex& x)
{
        if (is_exactly_a<numeric>(n_)
            and is_exactly_a<numeric>(x)
            and (ex_to<numeric>(n_).info(info_flags::inexact)
                 or ex_to<numeric>(x).info(info_flags::inexact)))
                return chebyshev_T(n_, x).evalf();

        if (x.is_one())
                return x;
        if (x.is_minus_one())
                return power(x, n_);
        if (x.is_zero())
                return cos(Pi*n_/_ex2);
        if (not is_exactly_a<numeric>(n_)
            or not n_.is_integer())
                return chebyshev_T(n_, x).hold();

        numeric numn = ex_to<numeric>(n_);
        if (numn.is_zero())
                return _ex1;
        if (numn.is_one())
                return x;
        if (numn.is_negative())
                numn = numn.negative();
        epvector vec;
        fmpz_poly_t p;
        fmpz_poly_init(p);
        fmpz_poly_chebyshev_t(p, static_cast<ulong>(numn.to_long()));
        int len = fmpz_poly_length(p);
        ex currx = _ex1;
        for (int i = 0; i<len; ++i) {
                mpz_t bigint;
                mpz_init(bigint);
                fmpz_get_mpz(bigint, fmpz_poly_get_coeff_ptr(p, i));
                numeric coeff(bigint);
                if (not coeff.is_zero())
                        vec.emplace_back(currx, coeff);
                currx *= x;
        }
        fmpz_poly_clear(p);
        return add(vec);
}

static ex chebyt_deriv(const ex& n, const ex & x, unsigned deriv_param)
{
	    if (deriv_param == 0)
                    throw std::runtime_error("derivative w.r.t. to the index is not supported yet");
	    return n * chebyshev_U(n-_ex1, x).hold();
}

REGISTER_FUNCTION(chebyshev_T, eval_func(chebyt_eval).
                        derivative_func(chebyt_deriv).
		        latex_name("T"));

//////////
// Chebyshev polynomials U_n(x)
//////////

static ex chebyu_eval(const ex& n_, const ex& x)
{
        if (is_exactly_a<numeric>(n_)
            and is_exactly_a<numeric>(x)
            and (ex_to<numeric>(n_).info(info_flags::inexact)
                 or ex_to<numeric>(x).info(info_flags::inexact)))
                return chebyshev_U(n_, x).evalf();

        if (x.is_one())
                return n_ + _ex1;
        if (x.is_minus_one())
                return power(_ex_1, n_) * (n_+_ex1);
        if (x.is_zero())
                return cos(Pi*n_/_ex2);
        if (not is_exactly_a<numeric>(n_)
            or not n_.is_integer())
                return chebyshev_U(n_, x).hold();

        const numeric& numn = ex_to<numeric>(n_);
        if (numn.is_zero())
                return _ex1;
        if (numn.is_one())
                return mul(x, _ex2);
        if (numn.is_negative())
                return -chebyu_eval(numn.negative()-*_num2_p, x);
        epvector vec;
        fmpz_poly_t p;
        fmpz_poly_init(p);
        fmpz_poly_chebyshev_u(p, static_cast<ulong>(numn.to_long()));
        int len = fmpz_poly_length(p);
        ex currx = _ex1;
        for (int i = 0; i<len; ++i) {
                mpz_t bigint;
                mpz_init(bigint);
                fmpz_get_mpz(bigint, fmpz_poly_get_coeff_ptr(p, i));
                numeric coeff(bigint);
                if (not coeff.is_zero())
                        vec.emplace_back(currx, coeff);
                currx *= x;
        }
        fmpz_poly_clear(p);
        return add(vec);
}

static ex chebyu_deriv(const ex& n, const ex & x, unsigned deriv_param)
{
	    if (deriv_param == 0)
                    throw std::runtime_error("derivative w.r.t. to the index is not supported yet");
	    return ((n+1) * chebyshev_T(n+_ex1, x).hold() - x*chebyshev_U(n, x)) / (power(x, 2) - _ex1);
}

REGISTER_FUNCTION(chebyshev_U, eval_func(chebyu_eval).
                        derivative_func(chebyu_deriv).
		        latex_name("U"));

//////////
// Legendre polynomials P_n(x)
//////////

static ex legp_evalf(const ex& n, const ex& x, PyObject* parent)
{
	if (not is_exactly_a<numeric>(x)
             or not is_exactly_a<numeric>(n))
                return legendre_P(n,x).hold();

        // see http://dlmf.nist.gov/15.9.E7
        const numeric& numn = ex_to<numeric>(n);
        const numeric& numx = ex_to<numeric>(x);
        std::vector<numeric> numveca, numvecb;
        numveca.push_back(numn.negative());
        numveca.push_back(numn + *_num1_p);
        numvecb.push_back(*_num1_p);
        return hypergeometric_pFq(numveca, numvecb, (*_num1_p - numx).mul(*_num1_2_p), parent);
}

static ex legp_eval(const ex& n_, const ex& x)
{
        ex n;
        if (n_.info(info_flags::negative))
                n = _ex_1 - n_;
        else
                n = n_;
	if (is_exactly_a<numeric>(x)) {
                const numeric& numx = ex_to<numeric>(x);
                if (numx.is_one())
                        return _ex1;
                if (numx.is_zero())
                        if (n.is_integer()) {
                                if (n.is_zero())
                                        return _ex1;
                                if (n.is_one())
                                        return x;
                                if (n.info(info_flags::odd))
                                        return _ex0;
                                if (n.info(info_flags::even)) {
                                        if (is_exactly_a<numeric>(n)) {
                                                const numeric& numn = ex_to<numeric>(n);
                                                return (numn + *_num_1_p).factorial()
                                                        / numn.mul(*_num1_2_p).factorial().pow_intexp(2)
                                                        * numn
                                                        / _num2_p->pow_intexp(numn.to_int())
                                                        * _num_1_p->pow_intexp(numn.mul(*_num1_2_p).to_int());
                                        }
                                        return gamma(n)
                                                / pow(gamma(n / _ex2), _ex2)
                                                / pow(_ex2, n - _ex2)
                                                / n
                                                * pow(_ex_1, n / _ex2);
                                }
                        }
                if (is_exactly_a<numeric>(n)
                    and (numx.info(info_flags::inexact)
                         or n.info(info_flags::inexact)))
                        return legp_evalf(n, x, nullptr);
        }

        if (not is_exactly_a<numeric>(n)
            or not n.is_integer())
                return legendre_P(n, x).hold();

        const numeric& numn = ex_to<numeric>(n);
        if (numn.is_zero())
                return _ex1;
        if (numn.is_one())
                return x;

        unsigned long d, k, L;
        unsigned long n_long = static_cast<unsigned long>(numn.to_long());
        d = k = n_long >> 1;

        while (k)
        {
                k >>= 1;
                d += k;
        }
        numeric den = _num2_p->pow_intexp(d);
        L = n_long / 2;
        unsigned int odd = n_long % 2;
        unsigned long index = odd;
        numeric curr_coeff = numeric::binomial(n_long, L);
        curr_coeff *= den;
        if (odd)
                curr_coeff *= L + 1;
        curr_coeff /= _num2_p->pow_intexp(2*L);
        if (L % 2)
                curr_coeff = curr_coeff.negative();

        epvector vec;
        vec.emplace_back(power(x, numeric(index)), curr_coeff / den);
        for (k = 1; k <= L; k++)
        {
                curr_coeff *= L + 1 - k;
                curr_coeff *= 2*k + 2*L - 1 + 2*odd;
                curr_coeff /= k;
                curr_coeff /= 2*k - 1 + 2*odd;
                curr_coeff = curr_coeff.negative();
                index += 2;
                vec.emplace_back(power(x, numeric(index)), curr_coeff / den);
        }

        return add(vec);
}

static ex legp_deriv(const ex& n, const ex & x, unsigned deriv_param)
{
	    if (deriv_param == 0)
                    throw std::runtime_error("derivative w.r.t. to the index is not supported yet");
	    return (n*legendre_P(n-1, x).hold() - n*x*legendre_P(n, x).hold()) / (1 - pow(x, 2));
}

REGISTER_FUNCTION(legendre_P, eval_func(legp_eval).
                        evalf_func(legp_evalf).
                        derivative_func(legp_deriv).
		        latex_name("P"));

//////////
// Hermite polynomials H_n(x)
//////////

static ex hermite_evalf(const ex& n, const ex& x, PyObject* parent)
{
	if (not is_exactly_a<numeric>(x)
             or not is_exactly_a<numeric>(n))
                return hermite(n,x).hold();

        // see http://dlmf.nist.gov/18.5.E13
        const numeric& numn = ex_to<numeric>(n);
        const numeric& numx = ex_to<numeric>(x);
        std::vector<numeric> numveca, numvecb;
        numveca.push_back(numn / *_num_2_p);
        numveca.push_back(*_num1_2_p + (numn / *_num_2_p));
        return pow(numx * (*_num2_p), numn) * hypergeometric_pFq(numveca, numvecb, numx.power(-2).negative(), parent);
}

static ex hermite_eval(const ex& n, const ex& x)
{
	if (is_exactly_a<numeric>(x)) {
                const numeric& numx = ex_to<numeric>(x);
                if (numx.is_zero())
                        if (n.is_integer() and n.info(info_flags::odd))
                                return _ex0;
                if (is_exactly_a<numeric>(n) and numx.info(info_flags::inexact))
                        return hermite_evalf(n, x, nullptr);
        }

        if (not is_exactly_a<numeric>(n))
                return hermite(n, x).hold();
        numeric numn = ex_to<numeric>(n);
        if (not numn.is_integer() or numn < 0)
                throw std::runtime_error("hermite_eval: The index n must be a nonnegative integer");

        if (numn.is_zero())
                return _ex1;

        // apply the formula from http://oeis.org/A060821
        // T(n, k) = ((-1)^((n-k)/2))*(2^k)*n!/(k!*((n-k)/2)!) if n-k is even and >=0, else 0.
        // sequentially for all viable k. Effectively there is the recurrence
        // T(n, k) = -(k+2)*(k+1)/(2*(n-k)) * T(n, k+2), with T(n, n) = 2^n
        numeric coeff = _num2_p->pow_intexp(numn);
        ex sum = _ex0;
        int fac = 1;
        epvector vec;
        while (numn >= 0) {
                vec.emplace_back(power(x, numn), coeff);
                coeff /= *_num_4_p;
                coeff *= numn;
                --numn;
                coeff *= numn;
                --numn;
                coeff /= fac++;
                }
        return add(vec);
}

static ex hermite_deriv(const ex& n, const ex & x, unsigned deriv_param)
{
	    if (deriv_param == 0)
                    throw std::runtime_error("derivative w.r.t. to the index is not supported yet");
	    return _ex2 * n * hermite(n-1, x).hold();
}

REGISTER_FUNCTION(hermite, eval_func(hermite_eval).
                        evalf_func(hermite_evalf).
                        derivative_func(hermite_deriv).
		        latex_name("H"));

//////////
// Gegenbauer (ultraspherical) polynomials C^a_n(x)
//////////

static ex gegenb_evalf(const ex& n, const ex &a, const ex& x, PyObject* parent)
{
	if (not is_exactly_a<numeric>(x)
             or not is_exactly_a<numeric>(a)
             or not is_exactly_a<numeric>(n))
                return gegenbauer(n,a,x).hold();

        // see http://dlmf.nist.gov/18.5.E9
        const numeric& numn = ex_to<numeric>(n);
        const numeric& numx = ex_to<numeric>(x);
        const numeric& numa = ex_to<numeric>(a);
        numeric num2a = numa * (*_num2_p);
        std::vector<numeric> numveca, numvecb;
        numveca.push_back(-numn);
        numveca.push_back(numn + num2a);
        numvecb.push_back(numa + (*_num1_2_p));
        numeric factor = numn + num2a - (*_num1_p);
        factor = factor.binomial(numn);
        return factor * hypergeometric_pFq(numveca, numvecb, ((*_num1_p)-numx)/(*_num2_p), parent);
}

static ex gegenb_eval(const ex& n, const ex &a, const ex& x)
{
	numeric numa = 1;

	if (is_exactly_a<numeric>(x)) {
                const numeric& numx = ex_to<numeric>(x);
                if (is_exactly_a<numeric>(n)
			and is_exactly_a<numeric>(a)
			and numx.info(info_flags::inexact))
                        return gegenb_evalf(n, a, x, nullptr);
        }

        if (not is_exactly_a<numeric>(n))
                return gegenbauer(n, a, x).hold();

        const numeric& numn = ex_to<numeric>(n);
        if (not numn.is_integer() or numn < 0)
                throw std::runtime_error("gegenb_eval: The index n must be a nonnegative integer");

        if (numn.is_zero())
                return _ex1;
	if (numn.is_equal(*_num1_p))
		return _ex2 * a * x;

	if (is_exactly_a<numeric>(a)) {
		numa = ex_to<numeric>(a);
		if (numa.is_zero())
			return _ex0; // C_n^0(x) = 0
	}

	if (not is_exactly_a<numeric>(a) or numa < 0) {
		ex p = _ex1;
		ex sum = _ex0;
		ex sign = _ex1;
		ex aa = _ex1;
		long nn = numn.to_long();

		sign = - 1;
		for (long k=0; k<=nn/2; k++) {
			sign *= - 1;
			aa = a;
			p = 1;

			// rising factorial (a)_{n-k}
			for (long i=0; i<nn-k; i++) {
				p *= aa;
				aa = aa + 1;
			}

			p /= numeric(k).factorial();
			p /= numeric(nn-2*k).factorial();
			sum += pow(2*x, nn-2*k) * p * sign;
		}
		return sum;
	}

	// from here on see flint2/fmpq_poly/gegenbauer_c.c
	numeric numer = numa.numer();
	numeric denom = numa.denom();
	numeric t = numn.factorial();
	numeric overall_denom = denom.pow_intexp(numn) * t;

	long nn = numn.to_long();
	numeric p = t / (numeric(nn/2).factorial());
	if ((nn%2) != 0u)
		p *= *_num2_p;
	if ((nn&2) != 0u)
		p = -p;

	for (long k=0; k < nn-nn/2; k++) {
		p *= numer;
		numer += denom;
	}

	p *= denom.power(nn/2);
        ex sum = _ex0;
	if ((nn%2) != 0u)
		sum += x*p;
	else
		sum += p;

	for (long k=nn/2-1; k>=0; --k) {
		p *= numer;
		p *= 4*(k+1);
		p *= *_num_1_p;
		p /= denom;
		p /= nn-2*k-1;
		p /= nn-2*k;
		sum += pow(x, nn-2*k) * p;
		numer += denom;
	}

        return sum / overall_denom;
}

static ex gegenb_deriv(const ex& n, const ex & a, const ex & x, unsigned deriv_param)
{
	    if (deriv_param == 0)
                    throw std::runtime_error("derivative w.r.t. to the index is not supported yet");
	    if (deriv_param == 1)
                    throw std::runtime_error("derivative w.r.t. to the second index is not supported yet");

	    return _ex2 * a * gegenbauer(n-1, a+1, x).hold();
}

REGISTER_FUNCTION(gegenbauer, eval_func(gegenb_eval).
                        evalf_func(gegenb_evalf).
                        derivative_func(gegenb_deriv).
		        latex_name("C"));

} // namespace GiNaC
