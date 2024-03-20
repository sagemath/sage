# sage_setup: distribution = sagemath-flint
cdef inline type element_class_zero = pAdicRelaxedElement_zero
cdef inline type element_class_one = pAdicRelaxedElement_one
cdef inline type element_class_bound = pAdicRelaxedElement_bound
cdef inline type element_class_value = pAdicRelaxedElement_value
cdef inline type element_class_random = pAdicRelaxedElement_random
cdef inline type element_class_slice = pAdicRelaxedElement_slice
cdef inline type element_class_add = pAdicRelaxedElement_add
cdef inline type element_class_sub = pAdicRelaxedElement_sub
cdef inline type element_class_mul = pAdicRelaxedElement_mul
cdef inline type element_class_muldigit = pAdicRelaxedElement_muldigit
cdef inline type element_class_div = pAdicRelaxedElement_div
cdef inline type element_class_sqrt = pAdicRelaxedElement_sqrt
cdef inline type element_class_teichmuller = pAdicRelaxedElement_teichmuller
cdef inline type element_class_unknown = pAdicRelaxedElement_unknown

include "sage/libs/linkages/padics/relaxed/flint.pxi"
include "relaxed_template.pxi"
