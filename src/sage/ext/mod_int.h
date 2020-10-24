/* See mod_int.pxd if you want to use these declarations in Cython */

/* CPython using signed long as the data type for plain integers, which is
 *   - 8 bytes on 64-bit Linux and OSX
 *   - 4 bytes on 64-bit Windows and all 32-bit supported platforms.
 * If you want to convert mod_int quickly to Python, then it must fit
 * into such an integer. Otherwise you'll end up with Python's
 * arbitrary-size integers instead of machine ints.
 *
 * However, 4 bytes is rather small and one can quickly run out of
 * primes, especially if the problem has many bad primes. See Trac
 * #10281 for details.
 *
 * Hence, we use signed 64-bit ints. This gives us fast conversion
 * to/from Python on 64-bit Linux and OSX, and a large number of
 * available (but not quite as fast) primes on 64-bit Windows and all
 * 32-bit platforms (see Trac #10281)
 */

typedef int64_t mod_int;

/* The largest value we can do arithmetic on without risking overflowing.
 * That is, you can multiply two MOD_INT_MAX in a mod_int.
 */
#define MOD_INT_OVERFLOW_UNSIGNED ((((uint64_t)1) << (sizeof(mod_int)*8 - 1)) - 1)
#define MOD_INT_OVERFLOW ((mod_int)MOD_INT_OVERFLOW_UNSIGNED)
#define MOD_INT_MAX ((mod_int)sqrt(MOD_INT_OVERFLOW_UNSIGNED) - 1)

