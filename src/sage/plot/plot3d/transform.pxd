# sage_setup: distribution = sagemath-plot
cdef struct point_c:
    double x, y, z

cdef struct color_c:
    double r, g, b

cdef struct face_c:
    int n
    int* vertices
    color_c color

cdef class Transformation:
    cdef matrix
    cdef double _matrix_data[12]
    cdef object _svd
    cpdef transform_point(self, x)
    cpdef transform_vector(self, v)
    cpdef transform_bounding_box(self, box)
    cdef void transform_point_c(self, point_c* res, point_c P) noexcept
    cdef void transform_vector_c(self, point_c* res, point_c P) noexcept
