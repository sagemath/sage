# sage_setup: distribution = sagemath-flint

# Deprecated header file; use sage/libs/flint/bernoulli.pxd instead
# See https://github.com/sagemath/sage/pull/36449

from sage.libs.flint.bernoulli cimport (
    bernoulli_fmpq_ui)
