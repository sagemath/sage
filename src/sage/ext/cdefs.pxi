#
# Declare C library functions used in Sage
#

from libc.stdio cimport *
from libc.string cimport strlen, strcpy, memset, memcpy, memcmp

from libc.math cimport sqrt, frexp, ldexp

from sage.libs.gmp.all cimport *
