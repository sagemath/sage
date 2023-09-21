/* sage_setup: distribution = sagemath-symbolics
 */
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
static struct sigaction ecl_sigfpe_handler;
static struct sigaction ecl_sigsegv_handler;
static struct sigaction sage_sigint_handler;
static struct sigaction sage_sigbus_handler;
static struct sigaction sage_sigfpe_handler;
static struct sigaction sage_sigsegv_handler;

static inline void set_ecl_signal_handler(void)
{
    sigaction(SIGINT, &ecl_sigint_handler, &sage_sigint_handler);
    sigaction(SIGBUS, &ecl_sigbus_handler, &sage_sigbus_handler);
    sigaction(SIGFPE, &ecl_sigfpe_handler, &sage_sigfpe_handler);
    sigaction(SIGSEGV, &ecl_sigsegv_handler, &sage_sigsegv_handler);
}

static inline void unset_ecl_signal_handler(void)
{
    sigaction(SIGINT, &sage_sigint_handler, NULL);
    sigaction(SIGBUS, &sage_sigbus_handler, NULL);
    sigaction(SIGFPE, &sage_sigfpe_handler, NULL);
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

#if ECL_VERSION_NUMBER < 230909
#define ecl_mpz_from_bignum(obj) ((obj)->big.big_num)
#else
#define ecl_mpz_from_bignum(obj) ecl_bignum(obj)
#endif

cl_object ecl_bignum_from_mpz(mpz_t num)
{
    cl_object z = _ecl_big_register0();
    mpz_set(ecl_mpz_from_bignum(z), num);
    return _ecl_big_register_copy(z);
}

static inline void safe_cl_boot(int argc, char** argv) {
   ECL_WITH_LISP_FPE_BEGIN {
      cl_boot(argc, argv);
   } ECL_WITH_LISP_FPE_END;
}

/* List of conditions to catch in the following functions. Is
 * initialized after cl_boot in init_ecl. */
static cl_object conditions_to_handle_clobj = ECL_NIL;

static inline cl_object safe_cl_funcall(cl_object *error, cl_object fun, cl_object arg) {
   cl_object ret = NULL;
   cl_env_ptr the_env = ecl_process_env();
   ECL_WITH_LISP_FPE_BEGIN {
       ECL_HANDLER_CASE_BEGIN(the_env, conditions_to_handle_clobj) {
           ret = cl_funcall(2, fun, arg);
       } ECL_HANDLER_CASE(1, condition) {
           *error = cl_princ_to_string(condition);
       } ECL_HANDLER_CASE_END;
   } ECL_WITH_LISP_FPE_END;
   return ret;
}

static inline cl_object safe_cl_apply(cl_object *error, cl_object fun, cl_object args) {
   cl_object ret = NULL;
   cl_env_ptr the_env = ecl_process_env();
   ECL_WITH_LISP_FPE_BEGIN {
       ECL_HANDLER_CASE_BEGIN(the_env, conditions_to_handle_clobj) {
           ret = cl_apply(2, fun, args);
       } ECL_HANDLER_CASE(1, condition) {
           *error = cl_princ_to_string(condition);
       } ECL_HANDLER_CASE_END;
   } ECL_WITH_LISP_FPE_END;
   return ret;
}

static inline cl_object safe_cl_eval(cl_object *error, cl_object form) {
   cl_object ret = NULL;
   cl_env_ptr the_env = ecl_process_env();
   ECL_WITH_LISP_FPE_BEGIN {
       ECL_HANDLER_CASE_BEGIN(the_env, conditions_to_handle_clobj) {
           ret = cl_eval(form);
       } ECL_HANDLER_CASE(1, condition) {
           *error = cl_princ_to_string(condition);
       } ECL_HANDLER_CASE_END;
   } ECL_WITH_LISP_FPE_END;
   return ret;
}
