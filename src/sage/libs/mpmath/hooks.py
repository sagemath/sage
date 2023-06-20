# sage_setup: distribution = sagemath-mpmath
r"""
Sage functions called by the vendored copy of mpmath

Upstream mpmath imports these functions from :mod:`sage.all`.

Our vendored copy of mpmath, :mod:`sage.libs.mpmath._vendor.mpmath`,
imports them from this module instead, see
``SAGE_ROOT/pkgs/sagemath-mpmath/pyproject.toml``.
"""
from sage.arith.misc import factorial, primes
from sage.combinat.combinat import fibonacci
from sage.rings.integer import Integer
