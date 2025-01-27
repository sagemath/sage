"""
Convert PARI objects to/from NTL objects

AUTHORS:

- Vincent Delecroix (2024)
"""
#*****************************************************************************
#       Copyright (C) 2024 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************


cdef Gen new_gen_from_ZZ(ZZ_c * x):
    sig_on()
    return new_gen(_new_GEN_from_ZZ(x))


