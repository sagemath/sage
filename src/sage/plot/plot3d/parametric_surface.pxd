from .index_face_set cimport IndexFaceSet
from .transform cimport point_c

cdef class ParametricSurface(IndexFaceSet):
    cdef object f
    cdef object render_grid
    cdef object color_function
    cdef object colormap
    cdef int eval_grid(self, urange, vrange) except -1
    cdef int eval_c(self, point_c *res, double u, double v) except -1

