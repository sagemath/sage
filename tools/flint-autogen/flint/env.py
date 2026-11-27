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


AUTOGEN_DIR = os.path.dirname(os.path.realpath(__file__))
FLINT_GIT_DIR = os.environ.get('FLINT_GIT_DIR', '')
FLINT_INCLUDE_DIR = os.path.join(FLINT_GIT_DIR, 'src')
FLINT_DOC_DIR = os.path.join(FLINT_GIT_DIR, 'doc/source')

if not os.path.isdir(FLINT_GIT_DIR) or not os.path.isdir(FLINT_INCLUDE_DIR) or not os.path.isdir(FLINT_DOC_DIR):
    raise ValueError('FLINT_GIT_DIR (={}) environment variable must be set to the location of flint sources'.format(FLINT_GIT_DIR))
