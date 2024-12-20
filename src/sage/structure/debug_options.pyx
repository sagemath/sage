# sage_setup: distribution = sagemath-objects
"""
Debug options for the :mod:`sage.structure` modules

EXAMPLES::

    sage: from sage.structure.debug_options import debug
    sage: debug.unique_parent_warnings
    False
    sage: debug.refine_category_hash_check
    True
    sage: debug.test_category_graph
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
        self.refine_category_hash_check = False
        self.test_category_graph = True
        self.test_nonrecursive_cachefunc = True

    def enable_extra_debugging_during_doctest(self):
        """
        Function that is called before doctest to enable extra debugging options.
        """
        self.refine_category_hash_check = True
        self.test_category_graph = True
        self.test_nonrecursive_cachefunc = True


cdef DebugOptions_class debug = DebugOptions_class()

# Since "debug" is declared with a type, it can only be cimported at
# the Cython level, not imported in plain Python. So we add it to the
# globals manually. We need to hack into Cython internals for this
# since Sage is compiled with the old_style_globals option.
from cpython.object cimport PyObject
cdef extern from *:
    PyObject* __pyx_d

(<object>__pyx_d)["debug"] = debug
