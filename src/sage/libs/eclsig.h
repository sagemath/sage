/*
 * eclsig.h: included from ecl.pyx
 *
 * AUTHORS:
 *
 * - Jeroen Demeyer (2011-03-30): initial version for #10818
 *
 */


#include <signal.h>
static struct sigaction ecl_sigint_handler;
static struct sigaction ecl_sigbus_handler;
static struct sigaction ecl_sigsegv_handler;
static struct sigaction sage_sigint_handler;
static struct sigaction sage_sigbus_handler;
static struct sigaction sage_sigsegv_handler;

static inline void set_ecl_signal_handler(void)
{
    sigaction(SIGINT, &ecl_sigint_handler, &sage_sigint_handler);
    sigaction(SIGBUS, &ecl_sigbus_handler, &sage_sigbus_handler);
    sigaction(SIGSEGV, &ecl_sigsegv_handler, &sage_sigsegv_handler);
}

static inline void unset_ecl_signal_handler(void)
{
    sigaction(SIGINT, &sage_sigint_handler, NULL);
    sigaction(SIGBUS, &sage_sigbus_handler, NULL);
    sigaction(SIGSEGV, &sage_sigsegv_handler, NULL);
}

/* This MUST be a macro because sig_on() must be in the same
 * stack frame as ecl_sig_on(). */
#define ecl_sig_on() \
    (sig_on() && (set_ecl_signal_handler() , 1))

static inline void ecl_sig_off(void)
{
    unset_ecl_signal_handler();
    sig_off();
}

#define ecl_mpz_from_bignum(obj) ((obj)->big.big_num)

cl_object ecl_bignum_from_mpz(mpz_t num)
{
    cl_object z = _ecl_big_register0();
    mpz_set(ecl_mpz_from_bignum(z), num);
    return _ecl_big_register_copy(z);
}
