# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_rng.h":
  cdef gsl_rng_type *gsl_rng_borosh13
  cdef gsl_rng_type *gsl_rng_coveyou
  cdef gsl_rng_type *gsl_rng_cmrg
  cdef gsl_rng_type *gsl_rng_fishman18
  cdef gsl_rng_type *gsl_rng_fishman20
  cdef gsl_rng_type *gsl_rng_fishman2x
  cdef gsl_rng_type *gsl_rng_gfsr4
  cdef gsl_rng_type *gsl_rng_knuthran
  cdef gsl_rng_type *gsl_rng_knuthran2
  cdef gsl_rng_type *gsl_rng_lecuyer21
  cdef gsl_rng_type *gsl_rng_minstd
  cdef gsl_rng_type *gsl_rng_mrg
  cdef gsl_rng_type *gsl_rng_mt19937
  cdef gsl_rng_type *gsl_rng_mt19937_1999
  cdef gsl_rng_type *gsl_rng_mt19937_1998
  cdef gsl_rng_type *gsl_rng_r250
  cdef gsl_rng_type *gsl_rng_ran0
  cdef gsl_rng_type *gsl_rng_ran1
  cdef gsl_rng_type *gsl_rng_ran2
  cdef gsl_rng_type *gsl_rng_ran3
  cdef gsl_rng_type *gsl_rng_rand
  cdef gsl_rng_type *gsl_rng_rand48
  cdef gsl_rng_type *gsl_rng_random128_bsd
  cdef gsl_rng_type *gsl_rng_random128_glibc2
  cdef gsl_rng_type *gsl_rng_random128_libc5
  cdef gsl_rng_type *gsl_rng_random256_bsd
  cdef gsl_rng_type *gsl_rng_random256_glibc2
  cdef gsl_rng_type *gsl_rng_random256_libc5
  cdef gsl_rng_type *gsl_rng_random32_bsd
  cdef gsl_rng_type *gsl_rng_random32_glibc2
  cdef gsl_rng_type *gsl_rng_random32_libc5
  cdef gsl_rng_type *gsl_rng_random64_bsd
  cdef gsl_rng_type *gsl_rng_random64_glibc2
  cdef gsl_rng_type *gsl_rng_random64_libc5
  cdef gsl_rng_type *gsl_rng_random8_bsd
  cdef gsl_rng_type *gsl_rng_random8_glibc2
  cdef gsl_rng_type *gsl_rng_random8_libc5
  cdef gsl_rng_type *gsl_rng_random_bsd
  cdef gsl_rng_type *gsl_rng_random_glibc2
  cdef gsl_rng_type *gsl_rng_random_libc5
  cdef gsl_rng_type *gsl_rng_randu
  cdef gsl_rng_type *gsl_rng_ranf
  cdef gsl_rng_type *gsl_rng_ranlux
  cdef gsl_rng_type *gsl_rng_ranlux389
  cdef gsl_rng_type *gsl_rng_ranlxd1
  cdef gsl_rng_type *gsl_rng_ranlxd2
  cdef gsl_rng_type *gsl_rng_ranlxs0
  cdef gsl_rng_type *gsl_rng_ranlxs1
  cdef gsl_rng_type *gsl_rng_ranlxs2
  cdef gsl_rng_type *gsl_rng_ranmar
  cdef gsl_rng_type *gsl_rng_slatec
  cdef gsl_rng_type *gsl_rng_taus
  cdef gsl_rng_type *gsl_rng_taus2
  cdef gsl_rng_type *gsl_rng_taus113
  cdef gsl_rng_type *gsl_rng_transputer
  cdef gsl_rng_type *gsl_rng_tt800
  cdef gsl_rng_type *gsl_rng_uni
  cdef gsl_rng_type *gsl_rng_uni32
  cdef gsl_rng_type *gsl_rng_vax
  cdef gsl_rng_type *gsl_rng_waterman14
  cdef gsl_rng_type *gsl_rng_zuf


  cdef gsl_rng_type *gsl_rng_default
  unsigned long int gsl_rng_default_seed

  gsl_rng *gsl_rng_alloc ( gsl_rng_type * T)
  int gsl_rng_memcpy (gsl_rng * dest, gsl_rng * src)
  gsl_rng *gsl_rng_clone ( gsl_rng * r)

  void gsl_rng_free (gsl_rng * r)

  void gsl_rng_set ( gsl_rng * r, unsigned long int seed)
  unsigned long int gsl_rng_max ( gsl_rng * r)
  unsigned long int gsl_rng_min ( gsl_rng * r)
  char *gsl_rng_name ( gsl_rng * r)

  int gsl_rng_fread (FILE * stream, gsl_rng * r)
  int gsl_rng_fwrite (FILE * stream, gsl_rng * r)

  size_t gsl_rng_size ( gsl_rng * r)
  void * gsl_rng_state ( gsl_rng * r)

  void gsl_rng_print_state ( gsl_rng * r)

  gsl_rng_type * gsl_rng_env_setup ()

  unsigned long int gsl_rng_get ( gsl_rng * r)
  double gsl_rng_uniform ( gsl_rng * r)
  double gsl_rng_uniform_pos ( gsl_rng * r)
  unsigned long int gsl_rng_uniform_int ( gsl_rng * r, unsigned long int n)


