# -*- coding: utf-8 -*-
"""
Environment Variables

This module defines the following subset of the Sage environment
variables:

* ``SAGE_ROOT``
* ``SAGE_SRC``
* ``SAGE_DISTFILES``
"""


# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os


try:
    SAGE_ROOT = os.environ['SAGE_ROOT']
except KeyError:
    SAGE_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))

SAGE_SRC = os.environ.get('SAGE_SRC',
                          os.path.join(SAGE_ROOT, 'src'))
SAGE_DISTFILES = os.environ.get('SAGE_DISTFILES',
                                os.path.join(SAGE_ROOT, 'upstream'))


assert os.path.isfile(os.path.join(SAGE_ROOT, 'configure.ac')), SAGE_ROOT

try:
    # SAGE_DISTFILES does not exist in a fresh git clone
    os.mkdir(SAGE_DISTFILES)
except OSError:
    pass

assert os.path.isdir(SAGE_DISTFILES)
