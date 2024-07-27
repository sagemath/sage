# sage_setup: distribution = sagemath-polyhedra
cimport cython
from sage.geometry.polyhedron.combinatorial_polyhedron.face_list_data_structure cimport face_list_t, face_t


@cython.final
cdef class ListOfFaces:
    # ``data`` points to the raw data.
    # It will be of "type" ``uint64_t[n_faces][face_length]``
    cdef face_list_t data

    cpdef ListOfFaces __copy__(self)

    cpdef int compute_dimension(self) except -2

    cdef inline size_t n_faces(self) noexcept:
        return self.data.n_faces
    cdef inline size_t n_atoms(self) noexcept:
        return self.data.n_atoms
    cdef inline size_t n_coatoms(self) noexcept:
        return self.data.n_coatoms

    cpdef ListOfFaces pyramid(self)

    cdef ListOfFaces delete_atoms_unsafe(self, bint* delete, face_t face)  # not in place
    cdef void delete_faces_unsafe(self, bint* delete, face_t face) noexcept  # in place

    cdef void get_not_inclusion_maximal_unsafe(self, bint *not_inclusion_maximal) noexcept
    cdef void get_faces_all_set_unsafe(self, bint *all_set) noexcept

cdef tuple face_as_combinatorial_polyhedron(ListOfFaces facets, ListOfFaces Vrep, face_t face, bint dual)
