cdef extern from "mag.h":
    ctypedef struct mag_struct:
        pass
    ctypedef mag_struct mag_t[1]
    ctypedef mag_struct * mag_ptr
    ctypedef const mag_struct * mag_srcptr
    long MAG_BITS

cdef extern from "arf.h":
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

cdef extern from "arb.h":
    ctypedef struct arb_struct:
        pass
    ctypedef arb_struct arb_t[1]
    ctypedef arb_struct * arb_ptr

cdef extern from "acb.h":
    ctypedef struct acb_struct:
        pass
    ctypedef acb_struct[1] acb_t
    ctypedef acb_struct * acb_ptr
    ctypedef const acb_struct * acb_srcptr

cdef extern from "acb_mat.h":
    ctypedef struct acb_mat_struct:
        pass
    ctypedef acb_mat_struct[1] acb_mat_t

cdef extern from "acb_poly.h":
    ctypedef struct acb_poly_struct:
        pass
    ctypedef acb_poly_struct[1] acb_poly_t
    ctypedef acb_poly_struct * acb_poly_ptr
    ctypedef const acb_poly_struct * acb_poly_srcptr
