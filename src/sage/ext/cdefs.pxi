"""
Deprecated include file

TESTS::

    sage: cython('include "sage/ext/cdefs.pxi"')
    doctest:...: DeprecationWarning: the file "cdefs.pxi" is deprecated, cimport the functions that you need
    See http://trac.sagemath.org/23855 for details.
"""

from sage.misc.superseded import deprecation
deprecation(23855, 'the file "cdefs.pxi" is deprecated, cimport the functions that you need')


from libc.stdio cimport *
from libc.string cimport strlen, strcpy, memset, memcpy, memcmp

from libc.math cimport sqrt, frexp, ldexp

from sage.libs.gmp.all cimport *
