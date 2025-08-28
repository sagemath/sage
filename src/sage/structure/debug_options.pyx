# sage_setup: distribution = sagemath-objects
"""
Debug options for the :mod:`sage.structure` modules

EXAMPLES::

    sage: from sage.structure.debug_options import debug
    sage: debug.unique_parent_warnings
    False
    sage: debug.refine_category_hash_check
    True
"""

#*****************************************************************************
#       Copyright (C) 2013 Robert Bradshaw <robertwb@gmail.com>
#                          William Stein <wstein@gmail.com>
#                          Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


cdef class DebugOptions_class:
    def __cinit__(self):
        """
        Initializer for the debug options.

        TESTS::

            sage: from sage.structure.debug_options import debug
            sage: type(debug)
            <... 'sage.structure.debug_options.DebugOptions_class'>
        """
        self.unique_parent_warnings = False
        # This one will be enabled during doctests
        self.refine_category_hash_check = False


cdef DebugOptions_class debug = DebugOptions_class()

# Since "debug" is declared with cdef, it can only be cimported at
# the Cython level, not imported in plain Python. So we add it to the
# globals manually. However this will make the variable out of sync
# if some user modifies the object, which is inevitable.
# See https://github.com/cython/cython/issues/3959#issuecomment-753455240
# and https://github.com/cython/cython/issues/656
# Note that ``_this_module`` could not be ``globals()``
# because Sage is compiled with the old_style_globals option.
import sage.structure.debug_options as _this_module
_this_module.debug = debug
del _this_module
