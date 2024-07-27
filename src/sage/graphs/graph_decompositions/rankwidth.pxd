# sage_setup: distribution = sagemath-graphs
cdef extern from "rw.h":
    ctypedef int uint_fast8_t
    ctypedef int uint_fast32_t
    ctypedef int subset_t
    int init_rw_dec(uint_fast8_t n)
    void destroy_rw()
    void calculate_level(uint_fast8_t subset_size)
    uint_fast8_t get_rw()
    subset_t *cslots
    subset_t *adjacency_matrix

cdef void print_rank_dec(subset_t s, int l) noexcept
