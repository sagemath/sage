#ifdef __cplusplus
/* The order here is very important. */
#include "eclib/curve.h"
#include "eclib/egr.h"
#include "eclib/descent.h"
#include "eclib/points.h"
#include "eclib/isogs.h"
#include "eclib/marith.h"
#endif

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

/**************** Miscellaneous functions ****************/

EXTERN long mwrank_get_precision(void);
EXTERN void mwrank_set_precision(long n);

EXTERN void mwrank_initprimes(char *pfilename, int verb);

/**************** bigint ****************/
#ifndef __cplusplus
struct bigint;
#endif

#ifdef __cplusplus
extern "C"
#endif
struct bigint* new_bigint(void);

EXTERN void del_bigint(struct bigint* x);

EXTERN struct bigint* str_to_bigint(char* s);


EXTERN char* bigint_to_str(struct bigint* x);


/**************** Curvedata ****************/

#ifndef __cplusplus
struct Curvedata;
#endif

EXTERN struct Curvedata* Curvedata_new(const struct bigint* a1, const struct bigint* a2,
                                       const struct bigint* a3, const struct bigint* a4,
                                       const struct bigint* a6, int min_on_init);

EXTERN void Curvedata_del(struct Curvedata* curve);


EXTERN char* Curvedata_repr(struct Curvedata* curve);

EXTERN double Curvedata_silverman_bound(const struct Curvedata* curve);

EXTERN double Curvedata_cps_bound(const struct Curvedata* curve);

EXTERN double Curvedata_height_constant(const struct Curvedata* curve);

EXTERN char* Curvedata_getdiscr(struct Curvedata* curve);

EXTERN char* Curvedata_conductor(struct Curvedata* curve);

EXTERN char* Curvedata_isogeny_class(struct Curvedata* E, int verbose);

/**************** mw -- subgroup of the mordell-weil group ****************/
#ifndef __cplusplus
struct mw;
#endif

EXTERN struct mw* mw_new(struct Curvedata* curve, int verb, int pp, int maxr);

EXTERN void mw_del(struct mw* m);

EXTERN int mw_process(struct Curvedata* curve, struct mw* m,
                      const struct bigint* x, const struct bigint* y,
                      const struct bigint* z, int sat);

EXTERN char* mw_getbasis(struct mw* m);

EXTERN char* mw_regulator(struct mw* m);

EXTERN int mw_rank(struct mw* m);

/* Returns index and unsat long array, which user must deallocate */
EXTERN int mw_saturate(struct mw* m, struct bigint* index, char** unsat,
                       long sat_bd, int odd_primes_only);

EXTERN void mw_search(struct mw* m, char* h_lim, int moduli_option, int verb);



/**************** two_descent ****************/
#ifndef __cplusplus
struct two_descent;
#endif

EXTERN struct two_descent* two_descent_new(struct Curvedata* curve,
                                    int verb, int sel,
                                    long firstlim, long secondlim,
                                    long n_aux, int second_descent);

EXTERN void two_descent_del(struct two_descent* t);

EXTERN long two_descent_get_rank(struct two_descent* t);

EXTERN long two_descent_get_rank_bound(struct two_descent* t);

EXTERN long two_descent_get_selmer_rank(struct two_descent* t);

EXTERN char* two_descent_get_basis(struct two_descent* t);

EXTERN int two_descent_ok(const struct two_descent* t);

EXTERN long two_descent_get_certain(const struct two_descent* t);

EXTERN void two_descent_saturate(struct two_descent* t, long sat_bd); // = -1 for default set in eclib/src/qcurves/saturate.h (currently 100)

EXTERN char* two_descent_regulator(struct two_descent* t);

