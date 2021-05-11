cimport cython
from sage.ext.memory_allocator  cimport MemoryAllocator
from sage.structure.sage_object cimport SageObject
from .list_of_faces             cimport ListOfFaces
from .face_data_structure       cimport face_t
from .face_iterator             cimport FaceIterator

@cython.final
cdef class CombinatorialFace(SageObject):
    cdef readonly bint _dual        # if 1, then iterate over dual Polyhedron
    cdef face_t face                # the face in bit-rep

    cdef MemoryAllocator _mem
    cdef size_t *atom_rep           # a place where atom-representation of face will be stored
    cdef size_t *coatom_rep         # a place where coatom-representation of face will be stored
    cdef int _dimension             # dimension of current face, dual dimension if ``dual``
    cdef int _ambient_dimension     # dimension of the polyhedron

    # An index to give different hashes for all faces of a Polyhedron.
    # The index must be chosen such that `F \subset G` implies ``hash(F) < hash(G)``.
    cdef size_t _hash_index

    cdef bint _initialized_from_face_lattice

    # some copies from ``CombinatorialPolyhedron``
    cdef tuple _ambient_Vrep, _ambient_facets, _equalities

    # Atoms and coatoms are the vertices/facets of the Polyedron.
    # If ``dual == 0``, then coatoms are facets, atoms vertices and vice versa.
    cdef ListOfFaces atoms, coatoms

    cdef size_t n_atom_rep(self) except -1
    cdef size_t set_coatom_rep(self) except -1
    cdef size_t set_atom_rep(self) except -1
