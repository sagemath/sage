from sage.plot.plot3d.parametric_surface cimport ParametricSurface


cdef class Cone(ParametricSurface):
    cdef double radius
    cdef double height
    cdef bint closed


cdef class Cylinder(ParametricSurface):
    cdef double radius
    cdef double height
    cdef bint closed


cdef class Sphere(ParametricSurface):
    cdef double radius


cdef class Torus(ParametricSurface):
    cdef double R, r
