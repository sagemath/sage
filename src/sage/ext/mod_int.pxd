# sage_setup: distribution = sagemath-objects
"""
The `mod_int` Data Type

* In C/C++ headers, you can `#include "mod_int.h"`
* In Cython files, use `from sage.ext.mod_int cimport *`
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2013 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


cdef extern from "mod_int.h":
    ctypedef long mod_int
    mod_int MOD_INT_MAX
    mod_int MOD_INT_OVERFLOW
