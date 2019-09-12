cimport cython
from libc.stdint                cimport uint64_t
from sage.ext.memory_allocator  cimport MemoryAllocator
from .list_of_faces             cimport ListOfFaces
from .combinatorial_face        cimport CombinatorialFace

@cython.final
cdef class PolyhedronFaceLattice:
    cdef MemoryAllocator _mem
    cdef int dimension              # dimension of Polyhedron
    cdef readonly bint dual         # if True, then List of all faces by dual Polyhedron
    cdef size_t face_length         # stores length of the faces in terms of uint64_t
    cdef tuple _V, _H, _equalities  # some copies from CombinatorialPolyhedron
    cdef size_t *f_vector           # a copy of the f-vector, is reversed if dual
    cdef size_t *face_counter       # how many faces of each dimension have been initialized
    cdef size_t *atom_repr          # a place where atom-representaion of face will be stored
    cdef size_t *coatom_repr        # a place where coatom-representaion of face will be stored

    # Atoms and coatoms are the Vrepr/facets of the Polyedron.
    # If ``dual == 0``, then coatoms are facets, atoms Vrepresentatives and vice versa.
    cdef ListOfFaces atoms, coatoms

    # All faces are stored in ``faces``. ``faces[i]`` stores all faces of
    # dimension `i` in Bit-representation (of atoms).
    cdef uint64_t *** faces
    cdef tuple faces_mem  # tuple to hold the ListOfFaces corresponding to faces

    # It follows a number of attributes to iterate over all incidences.
    # After initialization, ``PolyhedronFaceLattice`` can iterate over all incidences
    # of faces of dimension ``incidence_dim_one`` and ``incidence_dim_two``.
    cdef int is_incidence_initialized
    cdef int incidence_dim_one
    cdef int incidence_dim_two
    cdef size_t incidence_counter_one  # walks trough faces of incidence_dim_one
    cdef size_t incidence_counter_two  # walks through all indices of coatoms (for each face in incidence_dim_one)

    # Intersection of ``faces[incidence_dim_one][incidence_counter_one]`` with
    # ``coatoms[incidence_counter_two]``.
    # If this is contained in ``faces[incidence_dim_two]``, there is an incidence.
    cdef uint64_t *incidence_face

    cdef int _add_face(self, int face_dim, uint64_t *face) except -1
    cdef int _sort(self) except -1
    cdef int _sort_one_list(self, uint64_t **faces, size_t n_faces) except -1
    cdef int _sort_one_list_loop(
            self, uint64_t **inp, uint64_t **output1,
            uint64_t **output2, size_t n_faces) except -1
    cdef inline size_t find_face(self, int dimension, uint64_t *face) except -1
    cdef inline bint is_smaller(self, uint64_t *one, uint64_t *two)
    cdef inline int is_equal(self, int dimension, size_t index,
                             uint64_t *face) except -1
    cdef CombinatorialFace get_face(self, int dimension, size_t index)
    cdef size_t set_coatom_repr(self, int dimension, size_t index) except -1
    cdef size_t set_atom_repr(self, int dimension, size_t index) except -1
    cdef void incidence_init(self, int dimension_one, int dimension_two)
    cdef inline bint next_incidence(self, size_t *one, size_t *two)
    cdef inline bint next_incidence_loop(self, size_t *one, size_t *two)
    cdef inline bint next_trivial_incidence(self, size_t *one, size_t *two)
    cdef inline bint next_trivial_incidence2(self, size_t *one, size_t *two)
