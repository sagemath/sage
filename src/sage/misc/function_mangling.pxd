# sage_setup: distribution = sagemath-objects
cdef class ArgumentFixer:
    cdef public object f
    cdef public int _ndefault
    cdef public int _nargs
    cdef tuple _arg_names
    cdef bint _classmethod
    cdef dict _defaults
    cdef public tuple _default_tuple

    cdef fix_to_pos_args_kwds(self, tuple args, dict kwargs)
