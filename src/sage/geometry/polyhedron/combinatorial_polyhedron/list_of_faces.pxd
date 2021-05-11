cimport cython
from sage.ext.memory_allocator cimport MemoryAllocator
from .face_list_data_structure cimport face_list_t

@cython.final
cdef class ListOfFaces:
    cdef MemoryAllocator _mem

    # ``data`` points to the raw data.
    # It will be of "type" ``uint64_t[n_faces][face_length]``
    cdef face_list_t data

    cpdef int compute_dimension(self) except -2

    cdef inline size_t n_faces(self):
        return self.data.n_faces
    cdef inline size_t n_atoms(self):
        return self.data.n_atoms
    cdef inline size_t n_coatoms(self):
        return self.data.n_coatoms

    cpdef ListOfFaces pyramid(self)
