"""
Standard C helper code for Cython modules

Standard useful stuff for Sage Cython modules to include:
See stdsage.h for macros and stdsage.c for C functions.

Each module currently gets its own copy of this, which is why
we call the initialization code below.
"""
#*****************************************************************************
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef extern from "stdsage.h":
    ctypedef void PyObject

    # Global tuple -- useful optimization
    void init_global_empty_tuple()
    object PY_NEW(object t)
    object PY_NEW_SAME_TYPE(object t)

    void* PY_TYPE(object o)
    bint PY_TYPE_CHECK(object o, object t)
    bint PY_TYPE_CHECK_EXACT(object o, object t)

    object IS_INSTANCE(object o, object t)
    void PY_SET_TP_NEW(object t1, object t2)
    bint HAS_DICTIONARY(object o)
    bint PY_IS_NUMERIC(object o)


# Memory management
cdef extern from "stdsage.h":
    void  sage_free(void *p) nogil
    void* sage_realloc(void *p, size_t n) nogil
    void* sage_malloc(size_t) nogil
    void* sage_calloc(size_t nmemb, size_t size) nogil
    void  init_memory_functions() nogil
    void  init_csage()
    void  init_csage_module()


# Do this for every single module that links in stdsage.
init_csage_module()
