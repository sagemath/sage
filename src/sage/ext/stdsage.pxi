"""
Deprecated C helper code for Cython modules

TESTS::

    sage: cython('include "sage/ext/stdsage.pxi"')
    doctest:...: DeprecationWarning: the file "stdsage.pxi" is deprecated, cimport the functions that you need
    See http://trac.sagemath.org/23855 for details.
"""
#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.superseded import deprecation
deprecation(23855, 'the file "stdsage.pxi" is deprecated, cimport the functions that you need')


include "memory.pxi"

from cysignals.memory cimport sig_malloc as sage_malloc
from cysignals.memory cimport sig_realloc as sage_realloc
from cysignals.memory cimport sig_calloc as sage_calloc
from cysignals.memory cimport sig_free as sage_free
from cysignals.memory cimport (
        check_allocarray, check_reallocarray,
        check_malloc, check_realloc, check_calloc)

from sage.ext.stdsage cimport PY_NEW, HAS_DICTIONARY
from sage.ext.memory import init_memory_functions
