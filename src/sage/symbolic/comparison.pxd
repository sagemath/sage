from .ginac cimport *
from sage.symbolic.expression cimport Expression


cpdef int print_order(lhs, rhs) except -2
cdef int print_order_c(Expression lhs, Expression rhs)

cpdef print_sorted(expression_list)

# cpdef int math_order_c(Expression lhs, Expression rhs) except -2

