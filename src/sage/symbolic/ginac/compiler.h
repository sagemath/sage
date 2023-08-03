#ifndef GINAC_COMPILER_DEP_HH
#define GINAC_COMPILER_DEP_HH

#ifdef __GNUC__
#define unlikely(cond) __builtin_expect((cond), 0)
#define likely(cond) __builtin_expect((cond), 1)
#else
#define unlikely(cond) (cond)
#define likely(cond) (cond)
#endif

#endif /* GINAC_COMPILER_DEP_HH */
