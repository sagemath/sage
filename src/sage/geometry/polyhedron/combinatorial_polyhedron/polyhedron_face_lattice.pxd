cimport cython
from sage.ext.memory_allocator  cimport MemoryAllocator
from .list_of_faces             cimport ListOfFaces
from .face_data_structure       cimport face_t
from .face_list_data_structure  cimport face_list_t
from .combinatorial_face        cimport CombinatorialFace

@cython.final
cdef class PolyhedronFaceLattice:
    cdef MemoryAllocator _mem
    cdef int dimension              # dimension of Polyhedron
    cdef readonly bint dual         # if True, then List of all faces by dual Polyhedron
    cdef size_t *f_vector           # a copy of the f-vector, is reversed if dual
    cdef size_t *atom_rep           # a place where atom-representation of face will be stored
    cdef size_t *coatom_rep         # a place where coatom-representation of face will be stored

    # some copies from CombinatorialPolyhedron
    cdef tuple _Vrep, _facet_names, _equalities

    # Atoms and coatoms are the Vrep/facets of the Polyedron.
    # If ``dual == 0``, then coatoms are facets, atoms Vrepresentatives and vice versa.
    cdef ListOfFaces atoms, coatoms

    # All faces are stored in ``faces``. ``faces[i]`` stores all faces of
    # dimension `i` in Bit-representation (of atoms).
    cdef face_list_t *faces

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
    cdef face_t incidence_face

    cdef int _sort(self) except -1
    cdef inline size_t find_face(self, int dimension, face_t face) except -2
    cpdef CombinatorialFace get_face(self, int dimension, size_t index)
    cdef size_t set_coatom_rep(self, int dimension, size_t index) except -1
    cdef size_t set_atom_rep(self, int dimension, size_t index) except -1
    cdef void incidence_init(self, int dimension_one, int dimension_two)
    cdef inline bint next_incidence(self, size_t *one, size_t *two)
    cdef inline bint next_incidence_loop(self, size_t *one, size_t *two)
    cdef inline bint next_trivial_incidence(self, size_t *one, size_t *two)
    cdef inline bint next_trivial_incidence2(self, size_t *one, size_t *two)
