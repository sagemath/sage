# distutils: depends = mag.h arf.h arb.h acb.h acb_mat.h acb_poly.h acb_calc.h

# mag.h
cdef extern from "arb_wrap.h":
    ctypedef struct mag_struct:
        pass
    ctypedef mag_struct mag_t[1]
    ctypedef mag_struct * mag_ptr
    ctypedef const mag_struct * mag_srcptr
    long MAG_BITS

# arf.h
cdef extern from "arb_wrap.h":
    ctypedef struct arf_struct:
        pass
    ctypedef arf_struct arf_t[1]
    ctypedef arf_struct * arf_ptr
    ctypedef const arf_struct * arf_srcptr
    cdef enum arf_rnd_t:
        ARF_RND_DOWN
        ARF_RND_UP
        ARF_RND_FLOOR
        ARF_RND_CEIL
        ARF_RND_NEAR
    long ARF_PREC_EXACT

# arb.h
cdef extern from "arb_wrap.h":
    ctypedef struct arb_struct:
        pass
    ctypedef arb_struct arb_t[1]
    ctypedef arb_struct * arb_ptr
    ctypedef const arb_struct * arb_srcptr

# acb.h
cdef extern from "arb_wrap.h":
    ctypedef struct acb_struct:
        pass
    ctypedef acb_struct[1] acb_t
    ctypedef acb_struct * acb_ptr
    ctypedef const acb_struct * acb_srcptr

# acb_mat.h
cdef extern from "arb_wrap.h":
    ctypedef struct acb_mat_struct:
        pass
    ctypedef acb_mat_struct[1] acb_mat_t

# acb_poly.h
cdef extern from "arb_wrap.h":
    ctypedef struct acb_poly_struct:
        pass
    ctypedef acb_poly_struct[1] acb_poly_t
    ctypedef acb_poly_struct * acb_poly_ptr
    ctypedef const acb_poly_struct * acb_poly_srcptr

# acb_calc.h
cdef extern from "arb_wrap.h":
    ctypedef struct acb_calc_integrate_opt_struct:
        long deg_limit
        long eval_limit
        long depth_limit
        bint use_heap
        int verbose
    ctypedef acb_calc_integrate_opt_struct acb_calc_integrate_opt_t[1]
    ctypedef int (*acb_calc_func_t)(acb_ptr out,
            const acb_t inp, void * param, long order, long prec)

# arb_poly.h
cdef extern from "arb_wrap.h":
    ctypedef struct arb_poly_struct:
        pass
    ctypedef arb_poly_struct[1] arb_poly_t
    ctypedef arb_poly_struct * arb_poly_ptr
    ctypedef const arb_poly_struct * arb_poly_srcptr

