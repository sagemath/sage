# distutils: libraries = glpk z gmp

cdef extern from "glpk.h":
    int glp_init_env()
    int glp_free_env()
    int glp_term_out(int flag)
    void glp_term_hook(int (*func)(void *info, const char *s), void *info)
    int glp_open_tee(const char *fname)
    int glp_close_tee()
    void glp_error_hook(void (*func)(void *info), void *info)
    int glp_at_error()
    void glp_mem_limit(int limit)
    void glp_mem_usage(int *count, int *cpeak, size_t *total, size_t *tpeak)
