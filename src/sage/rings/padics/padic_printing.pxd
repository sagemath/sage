
include "sage/ext/cdefs.pxi"

from sage.structure.sage_object cimport SageObject
from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement

cdef class pAdicPrinter_class(SageObject):
    cdef object ring
    cdef int mode
    cdef bint pos
    cdef object ram_name
    cdef object unram_name
    cdef object var_name
    cdef object sep
    cdef object alphabet
    cdef PowComputer_class prime_pow
    cdef bint base
    cdef pAdicPrinter_class old
    cdef long max_ram_terms
    cdef long max_unram_terms
    cdef long max_terse_terms

    cdef base_p_list(self, value, bint pos)
    cdef _repr_gen(self, pAdicGenericElement elt, bint do_latex, bint pos, int mode, pname)
    cdef _repr_spec(self, pAdicGenericElement elt, bint do_latex, bint pos, int _mode, bint paren, pname)
    cdef _print_list_as_poly(self, L, bint do_latex, polyname, long expshift, bint increasing)
    cdef _truncate_list(self, L, max_terms, zero)
    cdef _var(self, x, exp, do_latex)
    cdef _dot_var(self, x, exp, do_latex)
    cdef _co_dot_var(self, co, x, exp, do_latex)
    cdef _plus_ellipsis(self, bint do_latex)
    cdef _ellipsis(self, bint do_latex)
    cdef _terse_frac(self, a, v, u, ram_name, bint do_latex)
    cdef _print_unram_term(self, L, bint do_latex, polyname, long max_unram_terms, long expshift, bint increasing)
    cdef _print_term_of_poly(self, s, coeff, bint do_latex, polyname, long exp)
