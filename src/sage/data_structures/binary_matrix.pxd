from sage.data_structures.bitset cimport bitset_t

cdef struct binary_matrix_s:
    Py_ssize_t n_cols
    Py_ssize_t n_rows
    bitset_t* rows
    
ctypedef binary_matrix_s[1] binary_matrix_t
