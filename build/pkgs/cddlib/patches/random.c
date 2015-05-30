/* Copyright (C) 1991, 1996 Free Software Foundation, Inc.

   ----------------------------------------------------------
   Random numbers that return the same sequence on all platforms.
   Implementation taken from the GNU C Library. The same copyright
   applies.
   ----------------------------------------------------------

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */


#include <stdint.h>
#include "random.h"


/* x**31 + x**3 + 1.  */
#define	TYPE_3		3
#define	BREAK_3		128
#define	DEG_3		31
#define	SEP_3		3


static int32_t randtbl[DEG_3 + 1] =
  {
    TYPE_3,
    
    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318,
  };


struct random_data
  {
    int32_t *fptr;              /* Front pointer.  */
    int32_t *rptr;              /* Rear pointer.  */
    int32_t *state;             /* Array of state values.  */
    int rand_type;              /* Type of random number generator.  */
    int rand_deg;               /* Degree of random number generator.  */
    int rand_sep;               /* Distance between front and rear.  */
    int32_t *end_ptr;           /* Pointer behind state table.  */
  };

static struct random_data rng_state =
  {
    .fptr = &randtbl[SEP_3 + 1],
    .rptr = &randtbl[1],
    .state = &randtbl[1],
    .rand_type = TYPE_3,
    .rand_deg = DEG_3,
    .rand_sep = SEP_3,
    .end_ptr = &randtbl[sizeof (randtbl) / sizeof (randtbl[0])]
};


void portable_srand(unsigned int seed)
{
  int type;
  int32_t *state;
  long int i;
  int32_t word;
  int32_t *dst;
  int kc;

  /* printf("portable_srand: seed = %u\n", seed); */

  type = rng_state.rand_type;

  state = rng_state.state;
  /* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
  if (seed == 0)
    seed = 1;
  state[0] = seed;

  dst = state;
  word = seed;
  kc = rng_state.rand_deg;
  for (i = 1; i < kc; ++i)
    {
      /* This does:
	   state[i] = (16807 * state[i - 1]) % 2147483647;
	 but avoids overflowing 31 bits.  */
      long int hi = word / 127773;
      long int lo = word % 127773;
      word = 16807 * lo - 2836 * hi;
      if (word < 0)
	word += 2147483647;
      *++dst = word;
    }

  rng_state.fptr = &state[rng_state.rand_sep];
  rng_state.rptr = &state[0];
  kc *= 10;
  while (--kc >= 0)
    portable_rand();
}


int portable_rand(void)
{
  int32_t *fptr = rng_state.fptr;
  int32_t *rptr = rng_state.rptr;
  int32_t *end_ptr = rng_state.end_ptr;
  int32_t result;

  result = *fptr += *rptr;
  /* Chucking least random bit.  */
  result = (result >> 1) & 0x7fffffff;
  ++fptr;
  if (fptr >= end_ptr)
    {
      fptr = rng_state.state;
      ++rptr;
    }
  else
    {
      ++rptr;
      if (rptr >= end_ptr)
	rptr = rng_state.state;
    }
  rng_state.fptr = fptr;
  rng_state.rptr = rptr;

  /* printf("portable_rand: output = %d\n", result); */

  return result;
}
