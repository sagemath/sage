# sage_setup: distribution = sagemath-plot
from sage.plot.plot3d.base cimport PrimitiveObject
from sage.plot.plot3d.transform cimport point_c, face_c, color_c


cdef class IndexFaceSet(PrimitiveObject):
    cdef bint enclosed
    cdef bint global_texture
    cdef Py_ssize_t vcount, fcount, icount
    cdef int realloc(self, Py_ssize_t vcount, Py_ssize_t fcount, Py_ssize_t icount) except -1
    # array of {x,y,z}
    cdef point_c* vs
    # pointers into face_indices marking the beginning of each face
    cdef face_c* _faces
    # array used as storage for _faces[i].vertices
    cdef int* face_indices


cdef class FaceIter:
    cdef Py_ssize_t i
    cdef IndexFaceSet set


cdef class EdgeIter:
    cdef Py_ssize_t i, j
    cdef object seen
    cdef IndexFaceSet set


cdef class VertexIter:
    cdef Py_ssize_t i
    cdef IndexFaceSet set
