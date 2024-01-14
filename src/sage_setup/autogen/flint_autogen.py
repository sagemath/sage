r"""
Autogeneration of flint headers.
"""
#*****************************************************************************
#       Copyright (C) 2023 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
from flint import AUTOGEN_DIR, write_flint_cython_headers


OUTPUT_DIR = os.path.realpath(os.path.join(AUTOGEN_DIR, os.path.pardir, os.path.pardir, os.path.pardir, 'sage', 'libs', 'flint'))
write_flint_cython_headers(OUTPUT_DIR, documentation=False)
