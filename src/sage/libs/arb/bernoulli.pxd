# sage_setup: distribution = sagemath-flint
# distutils: libraries = gmp flint
# distutils: depends = bernoulli.h

from sage.libs.flint.bernoulli cimport (
    bernoulli_fmpq_ui)
