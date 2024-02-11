from cpython.array cimport array

cdef void reset_swap(int n, int *c, int *o) noexcept
cdef int next_swap(int n, int *c, int *o) noexcept
cpdef bint next_perm(array l) noexcept
cpdef map_to_list(array l, tuple values, int n) noexcept
cpdef list left_action_same_n(list l, list r) noexcept
cpdef list right_action_same_n(list l, list r) noexcept
cpdef list left_action_product(list l, list r) noexcept
cpdef list right_action_product(list l, list r) noexcept

