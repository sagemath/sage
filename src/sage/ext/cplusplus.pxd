# sage_setup: distribution = sagemath-objects
#*****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef extern from "ccobject.h":
    # Print representation of any C++ object
    str ccrepr[T](const T& x)

    # Read a Python bytes/str into a C++ object
    int ccreadstr[T](T x, object b) except -1
