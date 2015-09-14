"""
Standard C helper code for Cython modules
"""
#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/interrupt.pxi"

from sage.ext.stdsage cimport PY_NEW, HAS_DICTIONARY
from sage.ext.memory cimport (
        sage_free, sage_realloc, sage_malloc, sage_calloc,
        check_allocarray, check_reallocarray,
        check_malloc, check_realloc, check_calloc)
from sage.ext.memory import init_memory_functions
