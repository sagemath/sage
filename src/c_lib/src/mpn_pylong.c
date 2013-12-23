/* mpn <-> pylong conversion and "pythonhash" for mpn
 *
 * Author:  Gonzalo Tornaría <tornaria@math.utexas.edu>
 * Date:    March 2006
 * License: GPL v2 or later
 *
 * the code to change the base to 2^PyLong_SHIFT is based on the function
 * mpn_get_str from GNU MP, but the new bugs are mine
 *
 * this is free software: if it breaks, you get to keep all the pieces
 */

#include "mpn_pylong.h"

/* This code assumes that PyLong_SHIFT < GMP_NUMB_BITS */
#if PyLong_SHIFT >= GMP_NUMB_BITS
#error "Python limb larger than GMP limb !!!"
#endif

/* Use these "portable" (I hope) sizebits functions
 * We could implement this in terms of count_leading_zeros from GMP,
 * but it is not exported !
 */
static const
unsigned char
__sizebits_tab[128] =
{
  0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7
};

#if GMP_LIMB_BITS > 64
#error "word size > 64 unsupported"
#endif

static inline
unsigned long
mpn_sizebits(mp_ptr up, mp_size_t un) {
  unsigned long cnt;
  mp_limb_t x;
  if (un==0) return 0;
  cnt = (un - 1) * GMP_NUMB_BITS;
  x = up[un - 1];
#if GMP_LIMB_BITS > 32
  if ((x >> 32) != 0) { x >>= 32; cnt += 32; }
#endif
#if GMP_LIMB_BITS > 16
  if ((x >> 16) != 0) { x >>= 16; cnt += 16; }
#endif
#if GMP_LIMB_BITS > 8
  if ((x >>  8) != 0) { x >>=  8; cnt += 8; }
#endif
  return cnt + ((x & 0x80) ? 8 : __sizebits_tab[x]);
}

static inline
unsigned long
pylong_sizebits(digit *digits, py_size_t size) {
  unsigned long cnt;
  digit x;
  if (size==0) return 0;
  cnt = (size - 1) * PyLong_SHIFT;
  x = digits[size - 1];
#if PyLong_SHIFT > 32
  if ((x >> 32) != 0) { x >>= 32; cnt += 32; }
#endif
#if PyLong_SHIFT > 16
  if ((x >> 16) != 0) { x >>= 16; cnt += 16; }
#endif
#if PyLong_SHIFT > 8
  if ((x >>  8) != 0) { x >>=  8; cnt += 8; }
#endif
  return cnt + ((x & 0x80) ? 8 : __sizebits_tab[x]);
}


/* mpn -> pylong conversion */

int
mpn_pylong_size (mp_ptr up, mp_size_t un)
{
  return (mpn_sizebits(up, un) + PyLong_SHIFT - 1) / PyLong_SHIFT;
}

/* this is based from GMP code in mpn/get_str.c */

/* Assume digits points to a chunk of size size
 * where size >= mpn_pylong_size(up, un)
 */
void
mpn_get_pylong (digit *digits, py_size_t size, mp_ptr up, mp_size_t un)
{
  mp_limb_t n1, n0;
  mp_size_t i;
  int bit_pos;
  /* point past the allocated chunk */
  digit * s = digits + size;

  /* input length 0 is special ! */
  if (un == 0) {
    while (size) digits[--size]=0;
    return;
  }

  i = un - 1;
  n1 = up[i];
  bit_pos = size * PyLong_SHIFT - i * GMP_NUMB_BITS;

  for (;;)
    {
      bit_pos -= PyLong_SHIFT;
      while (bit_pos >= 0)
        {
          *--s = (n1 >> bit_pos) & PyLong_MASK;
          bit_pos -= PyLong_SHIFT;
        }
      if (i == 0)
        break;
      n0 = (n1 << -bit_pos) & PyLong_MASK;
      n1 = up[--i];
      bit_pos += GMP_NUMB_BITS;
      *--s = n0 | (n1 >> bit_pos);
    }
}

/* pylong -> mpn conversion */

mp_size_t
mpn_size_from_pylong (digit *digits, py_size_t size)
{
  return (pylong_sizebits(digits, size) + GMP_NUMB_BITS - 1) / GMP_NUMB_BITS;
}

void
mpn_set_pylong (mp_ptr up, mp_size_t un, digit *digits, py_size_t size)
{
  mp_limb_t n1, d;
  mp_size_t i;
  int bit_pos;
  /* point past the allocated chunk */
  digit * s = digits + size;

  /* input length 0 is special ! */
  if (size == 0) {
    while (un) up[--un]=0;
    return;
  }

  i = un - 1;
  n1 = 0;
  bit_pos = size * PyLong_SHIFT - i * GMP_NUMB_BITS;

  for (;;)
    {
      bit_pos -= PyLong_SHIFT;
      while (bit_pos >= 0)
        {
          d = (mp_limb_t) *--s;
          n1 |= (d << bit_pos) & GMP_NUMB_MASK;
          bit_pos -= PyLong_SHIFT;
        }
      if (i == 0)
        break;
      d = (mp_limb_t) *--s;
      /* add some high bits of d; maybe none if bit_pos=-PyLong_SHIFT */
      up[i--] = n1 | (d & PyLong_MASK) >> -bit_pos;
      bit_pos += GMP_NUMB_BITS;
      n1 = (d << bit_pos) & GMP_NUMB_MASK;
    }
  up[0] = n1;
}


/************************************************************/

/* Hashing functions */

/* This is a bad hash...
 * If we decide to give up pylong compatibility, we should research to
 * find a decent (but fast) hash
 *
 * Some pointers to start:
 * <http://www.isthe.com/chongo/tech/comp/fnv/>
 * <http://www.azillionmonkeys.com/qed/hash.html>
 * <http://burtleburtle.net/bob/hash/doobs.html>
 */
/*
 * for an mpz, this number has to be multiplied by the sign
 * also remember to catch -1 and map it to -2 !
 */
long
mpn_pythonhash (mp_ptr up, mp_size_t un)
{
    /* Simply add all limbs */
    mp_limb_t h = 0;
    mp_limb_t h0;
    mp_size_t i;
    for (i = 0; i < un; i++)
    {
        h0 = h;
        h += up[i];
        /* Add 1 on overflow */
        if (h < h0) h++;
    }
    return h;
}

