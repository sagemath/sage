#*****************************************************************************
#       Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.perm_gps.partn_ref.data_structures cimport *


cdef class PythonPartitionStack:
    cdef PartitionStack *c_ps
