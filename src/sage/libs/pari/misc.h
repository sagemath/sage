/*****************************************
   PARI misc macros and functions
 *****************************************/
#ifndef SAGE_LIBS_PARI_MISC_H
#define SAGE_LIBS_PARI_MISC_H

#include <pari/pari.h>


int factorint_withproof_sage(GEN* ans, GEN x, GEN cutoff) {
  /*
  Factors and proves that the prime factors are really prime.
  If any aren't an ERROR condition (signal) is raised.

  INPUT:
     x -- a t_INT
     cutoff -- only check for primality of numbers at least this large.
  */

  GEN F = factorint(x, 0);
  *ans = F;

  long i, l;
  if (lg(F) == 1) return F; // x = 1
  F = gel(F,1); l = lg(F);
  for (i = 1; i < l; i++) {
    GEN p = gel(F,i);
    if (mpcmp(p, cutoff) > 0 && !isprime(p)) {
      char *c, *d;
      c = GENtostr(x);
      d = GENtostr(p);
      fprintf(stderr, "***\nPARI's factor(%s): Found composite pseudoprime %s (very rare and exciting -- PLEASE REPORT!!)\n***\n",
                   c, d);
      fprintf(stderr, "Do not worry, SAGE will further factor the number until each factor is proven prime.\n");
      free(c);
      free(d);
      return 1;
    }
  }
  return 0;
}

#endif  /* SAGE_LIBS_PARI_MISC_H */
