# sage_setup: distribution = sagemath-objects
"""
String <-> bytes encoding/decoding

TESTS:

Check that this can be used outside of Sage (see :issue:`25549`)::

    sage: cython(                                                                       # needs sage.misc.cython
    ....: '''
    ....: from sage.cpython.string cimport char_to_str
    ....: print(char_to_str("hello world!"))
    ....: ''')
    hello world!
"""

# ****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.bray@lri.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sys


# Provide this as a shortcut to calling sys.getfilesystemencoding(), which
# after interpreter initialization is constant.
FS_ENCODING = sys.getfilesystemencoding()

# Functions in this module are implemented in the .pxd file for inlining.
