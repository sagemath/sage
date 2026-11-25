#pragma once

#include <signal.h>
#include <setjmp.h>
#include <stdexcept>
#include <gmp.h>
#include <eclib/newforms.h>
#if defined(_ECLIB_MODULUS_H) && SCALAR_OPTION == 3
#include <flint/NTL-interface.h>
#endif

namespace sage_eclib {

// External installations through 20250627 use the homspace API without modulus.h.
#ifdef _ECLIB_MODULUS_H

class Homspace : public ::homspace {
public:
    const int cuspidal;

    Homspace(long level, int sign, int cusp, int verbose)
        : ::homspace(level, default_modulus<scalar>(), sign, verbose), cuspidal(cusp) {}

    long level() const { return N; }
    mat heckeop(long p, int dual, int verbose) const {
        return ::homspace::heckeop(p, cuspidal, dual, verbose);
    }
    smat s_heckeop(long p, int dual, int verbose) const {
        return ::homspace::s_heckeop(p, cuspidal, dual, verbose);
    }
};

inline newforms* make_newforms(long level, int verbose) {
    return new newforms(level, default_modulus<scalar>(), verbose);
}

#else

class Homspace : public ::homspace {
public:
    using ::homspace::homspace;
    long level() const { return modulus; }
};

inline newforms* make_newforms(long level, int verbose) {
    return new newforms(level, verbose);
}

#endif

#if defined(_ECLIB_MODULUS_H) && SCALAR_OPTION >= 3  // multiprecision entries

inline void scalar_to_mpz(mpz_ptr out, const scalar& value) {
    fmpz_t z;
#if SCALAR_OPTION == 4
    value.get_fmpz(z);  // INT::get_fmpz initializes its output
#else
    fmpz_init(z);
    fmpz_set_ZZ(z, value);
#endif
    fmpz_get_mpz(out, z);
    fmpz_clear(z);
}

inline scalar scalar_from_mpz(mpz_srcptr value) {
    fmpz_t z;
    fmpz_init(z);
    fmpz_set_mpz(z, value);
#if SCALAR_OPTION == 4
    scalar result(z);
#else
    scalar result;
    fmpz_get_ZZ(result, z);
#endif
    fmpz_clear(z);
    return result;
}

#else  // machine-integer entries

inline void scalar_to_mpz(mpz_ptr out, scalar value) {
    mpz_set_si(out, value);
}

inline scalar scalar_from_mpz(mpz_srcptr value) {
    if (!(sizeof(scalar) == sizeof(int) ? mpz_fits_sint_p(value) : mpz_fits_slong_p(value)))
        throw std::overflow_error("integer does not fit in an eclib scalar");
    return mpz_get_si(value);
}

#endif

}  // namespace sage_eclib
