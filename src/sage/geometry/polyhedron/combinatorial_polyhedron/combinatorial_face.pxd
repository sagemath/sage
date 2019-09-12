cimport cython
from libc.stdint                cimport uint64_t
from sage.ext.memory_allocator  cimport MemoryAllocator
from sage.structure.sage_object cimport SageObject
from .list_of_faces             cimport ListOfFaces
from .face_iterator             cimport FaceIterator

@cython.final
cdef class CombinatorialFace(SageObject):
    cdef readonly bint _dual        # if 1, then iterate over dual Polyhedron
    cdef ListOfFaces face_mem       # constructing face
    cdef uint64_t *face             # the face in bit-repr

    cdef MemoryAllocator _mem
    cdef size_t *atom_repr          # a place where atom-representaion of face will be stored
    cdef size_t *coatom_repr        # a place where coatom-representaion of face will be stored
    cdef int _dimension             # dimension of current face, dual dimension if ``dual``
    cdef int _ambient_dimension     # dimension of the polyhedron
    cdef size_t face_length         # stores length of the faces in terms of uint64_t
    cdef tuple _V, _H, _equalities  # some copies from ``CombinatorialPolyhedron``
    cdef size_t _hash_index         # an index to give different hashes for all faces of a Polyhedron

    # Atoms and coatoms are the vertices/facets of the Polyedron.
    # If ``dual == 0``, then coatoms are facets, atoms vertices and vice versa.
    cdef ListOfFaces atoms, coatoms

    cdef size_t length_atom_repr(self) except -1
    cdef size_t set_coatom_repr(self) except -1
    cdef size_t set_atom_repr(self) except -1
