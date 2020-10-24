cimport cython
from libc.stdint                cimport uint64_t
from sage.ext.memory_allocator  cimport MemoryAllocator
from sage.structure.sage_object cimport SageObject
from .list_of_faces             cimport ListOfFaces
from .combinatorial_face        cimport CombinatorialFace

cdef struct iter_struct:
    bint dual                  # if 1, then iterate over dual Polyhedron
    uint64_t *face             # the current face of the iterator
    size_t *atom_rep           # a place where atom-representaion of face will be stored
    size_t *coatom_rep         # a place where coatom-representaion of face will be stored
    int current_dimension      # dimension of current face, dual dimension if ``dual``
    int dimension              # dimension of the polyhedron
    int output_dimension       # only faces of this (dual?) dimension are considered
    int lowest_dimension       # don't consider faces below this (dual?) dimension
    size_t _index              # this counts the number of seen faces, useful for hasing the faces
    size_t face_length         # stores length of the faces in terms of uint64_t

    # ``visited_all`` points to faces, of which we have visited all faces already.
    # The number of faces in ``visited_all` might depend on the current dimension:
    #     Consider we visit the facets A,B of some face F.
    #     We will first visit all faces of A and then add A to visited_all.
    #     Then we visit all faces of B and add B to visited_all.
    #     Then we have visited F completely.
    #     Instead of having A and B in ``visited_all`` we will point to F.
    #     In this way, we will append ``visited_all`` in lower dimension, but
    #     will ignore those changes when going up in dimension again.
    #     This is why the number of faces in ``visited_all``depends on dimension.
    uint64_t **visited_all
    size_t *n_visited_all

    # ``maybe_newfaces`` is where all possible facets of a face are stored.
    # In dimension ``dim`` when visiting all faces of some face,
    # the intersections with other faces are stored in ``newfaces2[dim]``.
    uint64_t ***maybe_newfaces

    # ``newfaces`` will point to those faces in ``maybe_newfaces``
    # that are of codimension 1 and not already visited.
    uint64_t ***newfaces
    size_t *n_newfaces  # number of newfaces for each dimension

    # After having visited a face completely, we want to add it to ``visited_all``.
    # ``first_dim[i]`` will indicate, wether there is one more face in
    # ``newfaces[i]`` then ``n_newfaces[i]`` suggests
    # that has to be added to ``visited_all``.
    # If ``first_time[i] == False``, we still need to
    # add ``newfaces[i][n_newfaces[i]]`` to ``visited_all``.
    bint *first_time

    # The number of elements in newfaces[current_dimension],
    # that have not been visited yet.
    size_t yet_to_visit


cdef class FaceIterator_base(SageObject):
    cdef iter_struct structure
    cdef readonly bint dual         # if 1, then iterate over dual Polyhedron
    cdef MemoryAllocator _mem
    cdef tuple newfaces_lists       # tuple to hold the ListOfFaces corresponding to maybe_newfaces

    # some copies from ``CombinatorialPolyhedron``
    cdef tuple _Vrep, _facet_names, _equalities
    cdef bint _bounded

    # Atoms and coatoms are the vertices/facets of the Polyedron.
    # If ``dual == 0``, then coatoms are facets, atoms vertices and vice versa.
    cdef ListOfFaces atoms, coatoms

    cdef inline CombinatorialFace next_face(self)
    cdef inline int next_dimension(self) except -1
    cdef inline int next_face_loop(self) except -1
    cdef size_t n_atom_rep(self) except -1
    cdef size_t set_coatom_rep(self) except -1
    cdef size_t set_atom_rep(self) except -1
    cdef int ignore_subsets(self) except -1

@cython.final
cdef class FaceIterator(FaceIterator_base):
    pass

@cython.final
cdef class FaceIterator_geom(FaceIterator_base):
    cdef int _trivial_faces     # Whether to yield the trivial faces.
    cdef object _requested_dim  # Dimension requested on init.
    cdef readonly object P      # The original polyhedron.

# Nogil definitions of crucial functions.

cdef int next_dimension(iter_struct *structptr) nogil except -1
cdef int next_face_loop(iter_struct *structptr) nogil except -1
cdef size_t n_atom_rep(iter_struct *structptr) nogil except -1
